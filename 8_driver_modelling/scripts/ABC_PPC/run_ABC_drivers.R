
####################################################################################################
#
# "Round 1" simulations are the simulations used to perform ABC sampling.
# 
# At the start of "round 1", we generate a large sample of pairs (parameter vector, simulated data set)
# by sampling parameter vectors from the "prior" distribution,
# and for each parameter vector we run a simulation of the model;
# by the running script: sim_prior_drivers.R (which uses the R package rsimpop);
#
# This sample (the output from the script: sim_prior_drivers.R)
# is the input for the present script (collect_prior_sims_drivers.R);
#
# In fact, the script: sim_prior_drivers.R;
# is run many times within an "job array" on the computing cluster (farm5),
# so that may separate output files are generated;
# The present script (collect_prior_sims_drivers.R) collects and pools together all these output files;
#
# This large pooled sample of pairs (parameter vector, simulated data set)
# is then the input for the script: run_ABC_drivers.R;
#
# The present script (run_ABC_drivers.R) applies an ABC regression method (using the abc package)
# to this large pooled sample.
# The output from this script is a sample of observations (parameter vectors) from the posterior distribution.
# 
####################################################################################################

## Run on farm5
args = commandArgs(TRUE)
args

source_directory = toString(args[1])
tree_file = toString(args[2])
obs_stats_directory = toString(args[3])
obs_stats_file_name = toString(args[4])
sim_file = toString(args[5])
output_directory = toString(args[6])
posterior_sample_file_name = toString(args[7])
n_sims_accept = as.integer(args[8])
n_sims_max = as.integer(args[9])
patient_ID = toString(args[10])
stat_set_ID = toString(args[11])
n_stat_set=as.integer(args[12])

begin_index=13
end_index=begin_index+n_stat_set-1
stat_name_vect=args[begin_index:end_index]

cat( "0.0. run_ABC_drivers.R main: stat_name_vect = ", stat_name_vect, "\n", sep = "\t" )
cat( "0.0. run_ABC_drivers.R main: length( stat_name_vect ) = ", length( stat_name_vect ), "\n", sep = "\t" )
# readLines("stdin",n=1)


library(ape)
library(stringr)
library(tidyr)
library(seqinr)
library(phytools)
library(abc)
library(ggplot2)
library(devtools)
library(MCMCglmm)
library(phangorn)
library(spam)
library(INLA)
library(ggtree)
library(rsimpop)

#Set working directory
# my_working_directory = output_directory

#Set time points
# tree_path = paste0("/lustre/scratch117/casm/team154/em16/KX001/KX001_4/trees_vaf_cut/tree_",Run_ID,"_standard_rho01.tree")
tree_path = tree_file
#details_path = paste0()
stat_time=c(0.2,0.4,0.6,0.8)

##Initialise the environment
##Source functions
source(paste0(source_directory,"/functions/filters.parallel_functions.R"))
source(paste0(source_directory,"/functions/tree_functions.R"))
source(paste0(source_directory,"/functions/modelling_functions.R"))
source(paste0(source_directory,"/functions/targeted_analysis_functions.R"))
setwd("/lustre/scratch119/casm/team154pc/ms56/fetal_HSC/treemut"); source("treemut.R")
setwd(output_directory)

####################################################################################################
#
# Define functions needed for summary stats
#
get_summary_stats=function(tree,time_points) {
  nodeheights <- nodeHeights(tree)
  sumstats_each_timepoint=lapply(time_points, function(time_point) {
    clades_post_cutoff=tree$edge[,2][nodeheights[,1] < time_point & !nodeheights[,2] < time_point]
    clade_sizes=sapply(clades_post_cutoff,function(node) {length(getTips(tree,node))})
    n_singletons=sum(clade_sizes==1)
    clade_size_a=sort(clade_sizes, decreasing = T)[1]
    clade_size_b=sort(clade_sizes, decreasing = T)[2]
    clade_size_c=sort(clade_sizes, decreasing = T)[3]
    clade_size_d=sort(clade_sizes, decreasing = T)[4]
    clade_size_e=sort(clade_sizes, decreasing = T)[5]
    return(data.frame(n_singletons=n_singletons,clade_size_a=clade_size_a,clade_size_b=clade_size_b,clade_size_c=clade_size_c,clade_size_d=clade_size_d, clade_size_e=clade_size_e))
  })
  df=Reduce(cbind,sumstats_each_timepoint)
  return(df)
}

get_ltt = function(tree,time_points) {
  nodeheights <- nodeHeights(tree)
  ltt_tree = sapply(time_points, function(x) {
    sum(nodeheights[,1] < x & !nodeheights[,2] < x)
  })
  return(ltt_tree)
}

get_coalescences = function(ltt) {
  coals=sapply(2:length(ltt), function(i) {return(ltt[i]-ltt[i-1])})
  return(coals)
}

####################################################################################################
# This function maps the values of a parameter "theta" (restricted to range [a,b])
# onto the real line, using the "logit" function:
#
# z = logit( p ) = ln( p/(1-p) )
#

param_to_z <- function( theta, a, b, epsilon )
{
	
	p = ( theta - a ) / ( b - a )
	
	if( p <= 0 )
	{
		p = epsilon
		
	}else if( 1 <= p )
	{
		p = 1 - epsilon
		
	}
	
	z = log(p) - log(1-p)
	
	return( z )
	
}

####################################################################################################
# This function maps the values of a variable "z" (anywhere on the real line)
# onto the parameter range [a,b], using the "logistic" function:
#
# p = exp(z)/(1+exp(z)) = 1/(1+exp(-z))
#

z_to_param <- function( z, a, b, epsilon )
{
	
	p = 1/( 1 + exp(-z) )
	
	if( p <= 0 )
	{
		p = epsilon
		
	}else if( 1 <= p )
	{
		p = 1 - epsilon
		
	}
	
	theta = a + ( ( b - a ) * p )
	
	return( theta )
	
}

####################################################################################################
epsilon = 0.0001
####################################################################################################

cols_params= c("number_drivers_per_year","gamma_shape","gamma_rate")
cols_stats= c("mean_muts_sim", "n_singletons_1","clade_size_a_1","clade_size_b_1","clade_size_c_1","clade_size_d_1","clade_size_e_1","n_singletons_2","clade_size_a_2","clade_size_b_2","clade_size_c_2","clade_size_d_2","clade_size_e_2", "n_singletons_3", "clade_size_a_3","clade_size_b_3","clade_size_c_3","clade_size_d_3","clade_size_e_3","n_singletons_4", "clade_size_a_4","clade_size_b_4","clade_size_c_4","clade_size_d_4","clade_size_e_4", "ltt_1", "ltt_2", "ltt_3", "ltt_4", "coal_1", "coal_2", "coal_3")

sim_col_names_vect = c( cols_params, cols_stats )

####################################################################################################
#
# Read in observed data for individual "patient_ID";
#

#Get the data tree in shape and define the summary statistics 
system("mkdir pdfs")
tree.data=read.tree(tree_path) #Import the observed data tree
#load(details_path) #Import the details matrix
#tree.data.SNV = get_subset_tree(tree = tree.data, details = details, v.field="Mut_type", value = "SNV") #Create tree of only SNVs
mean_muts_data=mean(get_mut_burden(tree.data)) #Calculate the mean mutations per cell stat before make tree ultrametric
tree.data<-make.ultrametric.tree(tree.data) #Make ultrametric

# pdf(paste0("pdfs/",ID, "_obs_tree.pdf"), height =30)
pdf(paste0("pdfs/",patient_ID, "_obs_tree.pdf"), height =30)
plot(tree.data, show.tip.label = FALSE ); abline(v=stat_time,col="red")
dev.off()

#ltt_data=get_ltt(tree.data.ABC,time_points = time_points); names(ltt_data)=paste("ltt",1:length(ltt_data),sep="_")
#coalescences_data=get_coalescences(ltt_data); names(coalescences_data)=paste("coals",1:length(coalescences_data),sep="_")
#fractional_ltt_increase=coalscences_data/ltt_data[1:(length(ltt_data)-1)]; names(fractional_ltt_increase)=paste("frac_increase",1:length(fractional_ltt_increase),sep="_")
summary_stats <- get_summary_stats(tree.data,time_points=stat_time)
colnames(summary_stats) <- c("n_singletons_1","clade_size_a_1","clade_size_b_1","clade_size_c_1","clade_size_d_1","clade_size_e_1","n_singletons_2","clade_size_a_2","clade_size_b_2","clade_size_c_2","clade_size_d_2","clade_size_e_2", "n_singletons_3", "clade_size_a_3","clade_size_b_3","clade_size_c_3","clade_size_d_3","clade_size_e_3","n_singletons_4", "clade_size_a_4","clade_size_b_4","clade_size_c_4","clade_size_d_4","clade_size_e_4")
ltt_data <- get_ltt(tree.data,time_points=stat_time); names(ltt_data)=paste("ltt",1:length(ltt_data),sep="_")
coal_data <- get_coalescences(ltt_data);names(coal_data)=paste("coals",1:length(coal_data),sep="_")
stats_data = cbind(mean_muts_data,summary_stats,t(ltt_data),t(coal_data))

obs_stats_table = as.data.frame( matrix( NA, nrow = 1, ncol = length( cols_stats ) ) )
names( obs_stats_table ) <- cols_stats

obs_stats_table[ 1, cols_stats ] <- unlist(stats_data)

obs_stats_file = paste0( obs_stats_directory, "/", obs_stats_file_name )
write.table( obs_stats_table, file = obs_stats_file , append = FALSE , quote = FALSE , sep = "\t" , row.names = FALSE , col.names = TRUE )

####################################################################################################
#
# Read in simulations;
# Remove simulations which have missing data;
#

sim_table = read.table( file = sim_file, sep = "\t", header=TRUE, row.names = NULL, stringsAsFactors = FALSE )

# sumstats_to_include = c( "n_singletons_1","clade_size_a_1","clade_size_b_1","clade_size_c_1","clade_size_d_1","clade_size_e_1","n_singletons_2","clade_size_a_2","clade_size_b_2","clade_size_c_2","clade_size_d_2","clade_size_e_2", "n_singletons_3", "clade_size_a_3","clade_size_b_3","clade_size_c_3","clade_size_d_3","clade_size_e_3", "ltt_1", "ltt_2", "ltt_3", "coal_1", "coal_2" )

sumstats_to_include = stat_name_vect

abc_cols = c( cols_params, sumstats_to_include )
complete_rows_index_vect = complete.cases( sim_table[ , abc_cols ] )
sim_table = sim_table[ complete_rows_index_vect, ]

####################################################################################################

##NOW DO THE ABC

#
# Specify input to abc(...) function;
#

parameters.sim=sim_table[ , cols_params ]

parameters.sim[[ "number_drivers_per_year" ]] = sapply( parameters.sim[[ "number_drivers_per_year" ]], function(x){ param_to_z( x, a =1, b = 200, epsilon ) } )
parameters.sim[[ "gamma_shape" ]] = sapply( parameters.sim[[ "gamma_shape" ]], function(x){ param_to_z( x, a = 0.1, b =2, epsilon ) } )
parameters.sim[[ "gamma_rate" ]] = sapply( parameters.sim[[ "gamma_rate" ]], function(x){ param_to_z( x, a = 5, b = 120, epsilon ) } )

sumstat=sim_table[ , sumstats_to_include ]

target= unlist( obs_stats_table[ 1, sumstats_to_include ] )

if( nrow( sim_table ) > n_sims_max )
{
  sim_table = sim_table[ 1:n_sims_max, ]
}
n_sims = nrow( sim_table ) # n_sims = nrow( parameters.sim )
tol = n_sims_accept / n_sims

# abc_out = abc(target = target, param = parameters.sim, sumstat = sumstat, tol = tol, method = "rejection")
# abc_out = abc(target = target, param = parameters.sim, sumstat = sumstat, tol = tol, method = "loclinear")
abc_out = abc(target = target, param = parameters.sim, sumstat = sumstat, tol = tol, method = "ridge")

###############################################################################################
#
# Save ABC output to files
#

unadj_z_table = as.data.frame(abc_out$unadj.values)
adj_z_table = as.data.frame(abc_out$adj.values)

residuals_table = as.data.frame(abc_out$residuals)
weights_table = as.data.frame(abc_out$weights)
dist_table = as.data.frame(abc_out$dist)
# stats_accept_table = as.data.frame(abc_out$ss)

unadj_posterior_sample_table = as.data.frame(abc_out$unadj.values)
posterior_sample_table = as.data.frame(abc_out$adj.values)


posterior_sample_table[[ "number_drivers_per_year" ]] = sapply( posterior_sample_table[[ "number_drivers_per_year" ]], function(x){ z_to_param( x, a =1, b = 200, epsilon ) } )
posterior_sample_table[[ "gamma_shape" ]] = sapply( posterior_sample_table[[ "gamma_shape" ]], function(x){ z_to_param( x, a = 0.1, b =2, epsilon ) } )
posterior_sample_table[[ "gamma_rate" ]] = sapply( posterior_sample_table[[ "gamma_rate" ]], function(x){ z_to_param( x, a = 5, b = 120, epsilon ) } )

# posterior_sample_file = paste0( output_directory, "/", "posterior_sample.txt" )
posterior_sample_file = paste0( output_directory, "/", posterior_sample_file_name )
write.table( posterior_sample_table, file = posterior_sample_file , append = FALSE , quote = FALSE , sep = "\t" , row.names = FALSE , col.names = TRUE )


unadj_posterior_sample_table[[ "number_drivers_per_year" ]] = sapply( unadj_posterior_sample_table[[ "number_drivers_per_year" ]], function(x){ z_to_param( x, a =1, b = 200, epsilon ) } )
unadj_posterior_sample_table[[ "gamma_shape" ]] = sapply( unadj_posterior_sample_table[[ "gamma_shape" ]], function(x){ z_to_param( x, a = 0.1, b =2, epsilon ) } )
unadj_posterior_sample_table[[ "gamma_rate" ]] = sapply( unadj_posterior_sample_table[[ "gamma_rate" ]], function(x){ z_to_param( x, a = 5, b = 120, epsilon ) } )

unadj_posterior_sample_file = paste0( output_directory, "/", "unadj_posterior_sample", ".txt" )
write.table( unadj_posterior_sample_table, file = unadj_posterior_sample_file , append = FALSE , quote = FALSE , sep = "\t" , row.names = FALSE , col.names = TRUE )


unadj_z_file = paste0( output_directory, "/", "unadj_z", ".txt" )
write.table( unadj_z_table, file = unadj_z_file , append = FALSE , quote = FALSE , sep = "\t" , row.names = FALSE , col.names = TRUE )

adj_z_file = paste0( output_directory, "/", "adj_z", ".txt" )
write.table( adj_z_table, file = adj_z_file , append = FALSE , quote = FALSE , sep = "\t" , row.names = FALSE , col.names = TRUE )


residuals_file = paste0( output_directory, "/", "residuals", ".txt" )
write.table( residuals_table, file = residuals_file , append = FALSE , quote = FALSE , sep = "\t" , row.names = FALSE , col.names = TRUE )

weights_file = paste0( output_directory, "/", "weights", ".txt" )
write.table( weights_table, file = weights_file , append = FALSE , quote = FALSE , sep = "\t" , row.names = FALSE , col.names = TRUE )

dist_file = paste0( output_directory, "/", "dist", ".txt" )
write.table( dist_table, file = dist_file , append = FALSE , quote = FALSE , sep = "\t" , row.names = FALSE , col.names = TRUE )

###############################################################################################
#
# Plot distributions;
#
posterior_table = posterior_sample_table

prior_table = sim_table[ , cols_params ]

posterior_table$type="post"
prior_table$type="prior"
posterior_and_prior_table = rbind( posterior_table, prior_table )

# pdf(paste0("pdfs/",ID, "_number_drivers_per_year.pdf"))
pdf(paste0("pdfs/",patient_ID, "_number_drivers_per_year.pdf"))
# ggplot(data_and_prior, aes(x=number_drivers_per_year, col=factor(type))) + geom_density() + theme_classic()
ggplot(posterior_and_prior_table, aes(x=number_drivers_per_year, col=factor(type))) + geom_density() + theme_classic()
dev.off()

# pdf(paste0("pdfs/",ID, "_gamma_shape.pdf"))
pdf(paste0("pdfs/",patient_ID, "_gamma_shape.pdf"))
# ggplot(data_and_prior, aes(x=gamma_shape, col=factor(type))) + geom_density() + theme_classic()
ggplot(posterior_and_prior_table, aes(x=gamma_shape, col=factor(type))) + geom_density() + theme_classic()
dev.off()

# pdf(paste0("pdfs/",ID, "_gamma_rate.pdf"))
pdf(paste0("pdfs/",patient_ID, "_gamma_rate.pdf"))
# ggplot(data_and_prior, aes(x= gamma_rate, col=factor(type))) + geom_density() + theme_classic()
ggplot(posterior_and_prior_table, aes(x= gamma_rate, col=factor(type))) + geom_density() + theme_classic()
dev.off()

# pdf(paste0("pdfs/",ID, "_rate_mutation.pdf"))
pdf(paste0("pdfs/",patient_ID, "_rate_mutation.pdf"))
# ggplot(data, aes(x= gamma_rate, y = number_drivers_per_year)) + geom_density_2d() + theme_classic()
ggplot(posterior_and_prior_table, aes(x= gamma_rate, y = number_drivers_per_year)) + geom_density_2d() + theme_classic()
dev.off()

# pdf(paste0("pdfs/",ID, "_shape_mutation.pdf"))
pdf(paste0("pdfs/",patient_ID, "_shape_mutation.pdf"))
# ggplot(data, aes(x= gamma_shape, y = number_drivers_per_year)) + geom_density_2d() + theme_classic()
ggplot(posterior_and_prior_table, aes(x= gamma_shape, y = number_drivers_per_year)) + geom_density_2d() + theme_classic()
dev.off()

# pdf(paste0("pdfs/",ID, "_rate_shape.pdf"))
pdf(paste0("pdfs/",patient_ID, "_rate_shape.pdf"))
# ggplot(data, aes(x= gamma_shape, y = gamma_rate)) + geom_density_2d() + theme_classic() 
ggplot(posterior_and_prior_table, aes(x= gamma_shape, y = gamma_rate)) + geom_density_2d() + theme_classic() 
dev.off()


#Calculate and print quantile results
# number_drivers_per_year_result_95CI=quantile(rej.KX001$unadj.values[,"number_drivers_per_year"],probs=c(0.025,0.975))
# number_drivers_per_year_result_80CI=quantile(rej.KX001$unadj.values[,"number_drivers_per_year"],probs=c(0.1,0.9))
number_drivers_per_year_result_95CI=quantile(posterior_table[,"number_drivers_per_year"],probs=c(0.025,0.975))
number_drivers_per_year_result_80CI=quantile(posterior_table[,"number_drivers_per_year"],probs=c(0.1,0.9))
number_drivers_per_year_result_95CI
number_drivers_per_year_result_80CI

# gamma_shape_result_95CI=quantile(rej.KX001$unadj.values[,"gamma_shape"],probs=c(0.025,0.975))
# gamma_shape_result_80CI=quantile(rej.KX001$unadj.values[,"gamma_shape"],probs=c(0.1,0.9))
gamma_shape_result_95CI=quantile(posterior_table[,"gamma_shape"],probs=c(0.025,0.975))
gamma_shape_result_80CI=quantile(posterior_table[,"gamma_shape"],probs=c(0.1,0.9))
gamma_shape_result_95CI
gamma_shape_result_80CI

# gamma_rate_result_95CI=quantile(rej.KX001$unadj.values[,"gamma_rate"],probs=c(0.025,0.975))
# gamma_rate_result_80CI=quantile(rej.KX001$unadj.values[,"gamma_rate"],probs=c(0.1,0.9))
gamma_rate_result_95CI=quantile(posterior_table[,"gamma_rate"],probs=c(0.025,0.975))
gamma_rate_result_80CI=quantile(posterior_table[,"gamma_rate"],probs=c(0.1,0.9))
gamma_rate_result_95CI
gamma_rate_result_80CI


posterior_and_prior_file = paste0( output_directory, "/", patient_ID, "_data_and_prior" )
write.table( posterior_and_prior_table, file = posterior_and_prior_file , append = FALSE , quote = FALSE , sep = "\t" , row.names = FALSE , col.names = TRUE )

###############################################################################################

cat( "done!", "\n" )

