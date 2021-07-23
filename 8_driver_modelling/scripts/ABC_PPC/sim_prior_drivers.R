
####################################################################################################
#
# "Round 1" simulations are the simulations used to perform ABC sampling.
# 
# At the start of "round 1", we generate a large sample of pairs (parameter vector, simulated data set)
# by sampling parameter vectors from the "prior" distribution,
# and for each parameter vector we run a simulation of the model;
# by the running the present script (sim_prior_drivers.R), which uses the R package rsimpop;
#
# In fact, the script: sim_prior_drivers.R;
# is run many times within an "job array" on the computing cluster (farm5),
# so that may separate output files are generated;
#
# These separate output files are the input for the script: collect_prior_sims_drivers.R;
# which collects and pools together all these output files;
#
####################################################################################################

## Run on farm5
args = commandArgs(TRUE)
args

##SET FILEPATHS
# RUNID = as.character(args[1])

source_directory = toString(args[1])
output_directory = toString(args[2])
prior_sample_file = toString(args[3])
use_prior_sample=as.logical(as.integer(args[4]))
donor_age=as.numeric(args[5])
donor_n_genomes=as.integer(args[6])
n_sims_per_job=as.integer(args[7])
sim_job_index=as.integer(args[8])
# rand_seed=as.integer(args[9])

cat( "0.0. sim_prior_drivers.R main: use_prior_sample = ", use_prior_sample, "\n", sep = "\t" )
cat( "0.0. sim_prior_drivers.R main: donor_age = ", donor_age, "\n", sep = "\t" )
cat( "0.0. sim_prior_drivers.R main: donor_n_genomes = ", donor_n_genomes, "\n", sep = "\t" )
# readLines("stdin",n=1)

library(ape)
library(stringr)
library(tidyr)
library(seqinr)
library(phytools)
library(ggplot2)
library(devtools)
library(MCMCglmm)
library(phangorn)
library(spam)
library(INLA)
library(ggtree)
library(rsimpop)
library(truncdist)


# Create the output directory:
if( !file.exists( output_directory ) )
{
	dir.create( output_directory )
	
}

#Set time points
stat_time=c(0.2,0.4,0.6,0.8)
run <- round(runif( n=1, min= 1e6, max = 2e6), digits = 0)
# my_working_directory = ("/lustre/scratch117/casm/team154/em16/simulations/rsimpop_drivers/rsimpop/")
# my_working_directory = output_directory

##Initialise the environment
##Source functions
source(paste0(source_directory,"/functions/filters.parallel_functions.R"))
source(paste0(source_directory,"/functions/tree_functions.R"))
source(paste0(source_directory,"/functions/modelling_functions.R"))
source(paste0(source_directory,"/functions/targeted_analysis_functions.R"))
setwd("/lustre/scratch119/casm/team154pc/ms56/fetal_HSC/treemut"); source("treemut.R");
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
#
# Create vectors of column names;
# Create data frame "sim_table";
# Sample from prior;
#

# nsim_per_job = 20
nsim_per_job = n_sims_per_job

fitness_threshold = 0.05
number_params=3
cols_params= c("number_drivers_per_year","gamma_shape","gamma_rate")
number_stats=32
cols_stats= c("mean_muts_sim", "n_singletons_1","clade_size_a_1","clade_size_b_1","clade_size_c_1","clade_size_d_1","clade_size_e_1","n_singletons_2","clade_size_a_2","clade_size_b_2","clade_size_c_2","clade_size_d_2","clade_size_e_2", "n_singletons_3", "clade_size_a_3","clade_size_b_3","clade_size_c_3","clade_size_d_3","clade_size_e_3","n_singletons_4", "clade_size_a_4","clade_size_b_4","clade_size_c_4","clade_size_d_4","clade_size_e_4", "ltt_1", "ltt_2", "ltt_3", "ltt_4", "coal_1", "coal_2", "coal_3")


sim_col_names_vect = c( cols_params, cols_stats )
sim_table = as.data.frame( matrix( NA, nrow = nsim_per_job, ncol = length( sim_col_names_vect ) ) )
names( sim_table ) <- sim_col_names_vect

if( isTRUE( use_prior_sample ) )
{
  
cat( "1.1. sim_prior_drivers.R main: use_prior_sample = ", use_prior_sample, "\n", sep = "\t" )
cat( "1.1. sim_prior_drivers.R main: isTRUE( use_prior_sample ) = ", isTRUE( use_prior_sample ), "\n", sep = "\t" )
  
  prior_table = read.table( file = prior_sample_file, sep = "\t", header=TRUE, stringsAsFactors = FALSE )
  
  row_index_vect = sample( x=nrow( prior_table ), size=nsim_per_job, replace = TRUE )
  
  sim_table[ , cols_params ] = prior_table[ row_index_vect, ]
  
}else
{
  
cat( "1.0. sim_prior_drivers.R main: use_prior_sample = ", use_prior_sample, "\n", sep = "\t" )
cat( "1.0. sim_prior_drivers.R main: isTRUE( use_prior_sample ) = ", isTRUE( use_prior_sample ), "\n", sep = "\t" )
  
  # sim_table[[ "number_drivers_per_year" ]] = round(runif( n=nsim_per_job, min =1, max = 120), digits = 0)
  sim_table[[ "number_drivers_per_year" ]] = round(runif( n=nsim_per_job, min =1, max = 200), digits = 0)
  sim_table[[ "gamma_shape" ]] = round(runif( n=nsim_per_job, min = 0.1, max =2), digits = 2)
  sim_table[[ "gamma_rate" ]] = round(runif( n=nsim_per_job, min = 5, max = 120), digits = 2)
  
}

####################################################################################################
#
# Run simulations;
#

sim_outfile = paste( output_directory, "/", "sim_", "job_index_", sim_job_index, "_run_", run, ".txt", sep = "" )

for(i in 1:nsim_per_job) {
  set.seed(run+i)
  
cat( "2.0. sim_prior_drivers.R main: i = ", i, "\n", sep = "\t" )
  
  number_drivers_per_year = sim_table[[ "number_drivers_per_year" ]][ i ]
  gamma_shape = sim_table[[ "gamma_shape" ]][ i ]
  gamma_rate = sim_table[[ "gamma_rate" ]][ i ]
  
  # gamma_shape = round(runif( n=1, min = 0.1, max =2), digits = 2)
  # gamma_rate = round(runif( n=1, min = 5, max = 120), digits = 2)
  # number_drivers_per_year = round(runif( n=1, min =1, max = 120), digits = 0)
  param_ID = paste(number_drivers_per_year, gamma_shape, gamma_rate, sep = "_")
  parameters = list(number_drivers_per_year = number_drivers_per_year, gamma_shape = gamma_shape, gamma_rate = gamma_rate)
  print(parameters)
  
  ##Function to generate gamma distribution based fitness
  genGammaFitness=function(fitness_threshold,shape,rate){
    function() rtrunc(n=1,a=fitness_threshold, b=Inf,"gamma",shape=shape,rate=rate)
  }
  fitnessGammaFn=genGammaFitness(fitness_threshold=fitness_threshold,shape = gamma_shape, rate=gamma_rate)
  
  ####################################################################################################
  
  ### Run the simulation for donor
  # dps = run_driver_process_sim( 0.1, 0.5/(365), target_pop_size = 1e5, nyears = donor_age, fitness= fitnessGammaFn, drivers_per_year = number_drivers_per_year )
  dps = run_driver_process_sim_KEVIN( 0.1, 0.5/(365), target_pop_size = 1e5, nyears = donor_age, fitness= fitnessGammaFn, drivers_per_year = number_drivers_per_year )
  
  ### Continue the simulation
  # dps = continue_driver_process_sim( dps, nyears = donor_age, fitness= fitnessGammaFn )
  # dps = continue_driver_process_sim_KEVIN( dps, nyears = donor_age, fitness= fitnessGammaFn )
  
  ####################################################################################################
  
  ###Get subsampled tree
  # dpst_KX001=get_subsampled_tree(dps_KX001,KX001_WGS_sequenced)
  dpst = get_subsampled_tree( dps, donor_n_genomes )
  
  # tree=get_elapsed_time_tree(dpst_KX001,mutrateperdivision = 1,backgroundrate = 16/365)
  tree = get_elapsed_time_tree( dpst, mutrateperdivision = 1, backgroundrate = 16/365 )
  
  mean_muts_sim=mean(get_mut_burden(tree))
  
  tree <- make.ultrametric.tree(tree)
  
  ####################################################################################################
  
  
  summary_stats <- get_summary_stats(tree,time_points=stat_time)
  colnames(summary_stats) <- c("n_singletons_1","clade_size_a_1","clade_size_b_1","clade_size_c_1","clade_size_d_1","clade_size_e_1","n_singletons_2","clade_size_a_2","clade_size_b_2","clade_size_c_2","clade_size_d_2","clade_size_e_2", "n_singletons_3", "clade_size_a_3","clade_size_b_3","clade_size_c_3","clade_size_d_3","clade_size_e_3","n_singletons_4", "clade_size_a_4","clade_size_b_4","clade_size_c_4","clade_size_d_4","clade_size_e_4")
  ltt_sim <- get_ltt(tree,time_points=stat_time); names(ltt_sim)=paste("ltt",1:length(ltt_sim),sep="_")
  coal_sim <- get_coalescences(ltt_sim);names(coal_sim)=paste("coal",1:length(coal_sim),sep="_")
  stats_sim = cbind(mean_muts_sim,summary_stats,t(ltt_sim),t(coal_sim))
  
  sim_table[ i, cols_stats ] <- unlist(stats_sim)
  
  if( i == 1 )
  {
	  write.table( sim_table[ i, ], file = sim_outfile , append = FALSE , quote = FALSE , sep = "\t" , row.names = FALSE , col.names = TRUE )
	  # write.table( sim_table[ i, ], file = sim_outfile , append = TRUE , quote = FALSE , sep = "\t" , row.names = FALSE , col.names = FALSE )
	  
  }else
  {
	  # write.table( sim_table[ i, ], file = sim_outfile , append = FALSE , quote = FALSE , sep = "\t" , row.names = FALSE , col.names = TRUE )
	  write.table( sim_table[ i, ], file = sim_outfile , append = TRUE , quote = FALSE , sep = "\t" , row.names = FALSE , col.names = FALSE )
	  
  }
  
} # for(i in 1:nsim_per_job)

###############################################################################################

cat( "done!", "\n" )
















