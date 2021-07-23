
####################################################################################################
#
# "Round 1" simulations are the simulations used to perform ABC sampling.
# "Round 2" simulations ("resimulations") are the simulations used to perform PPC calculations.
#
# After completion of ABC computation at the end of "round 1", we have
# a sample of parameter vectors from the (approximate) posterior distribution;
# This is the output from the script: run_ABC_drivers.R;
# 
# At the start of "round 2", we take this pre-existing sample of parameter vectors from this (approximate) posterior distribution
# (OR a sub-sample drawn withOUT replacment from this pre-existing sample of parameter vectors),
# and for each parameter vector in this sample, we perform (a large number of) simulations
# from conditional distribution at each parameter vector;
# by running the present script (resim_conditionals_drivers.R);
#
####################################################################################################

## Run on farm5
args = commandArgs(TRUE)
args

##SET FILEPATHS
# RUNID = as.character(args[1])

source_directory = toString(args[1])
output_path = toString(args[2])
posterior_sample_file = toString(args[3])
resim_output_directory_name = toString(args[4])
n_sims_per_job=as.integer(args[5])
post_obs_index=as.integer(args[6])
resim_job_index=as.integer(args[7])
n_donors=as.integer(args[8])

begin_index=9
end_index=begin_index+n_donors-1
donor_ID_vect=args[begin_index:end_index]

begin_index=begin_index+n_donors
end_index=begin_index+n_donors-1
donor_age_vect=args[begin_index:end_index]
donor_age_vect=as.integer(donor_age_vect)

begin_index=begin_index+n_donors
end_index=begin_index+n_donors-1
donor_n_genomes_vect=args[begin_index:end_index]
donor_n_genomes_vect=as.integer(donor_n_genomes_vect)

cat( "0.0. resim_conditionals_drivers.R main: donor_ID_vect = ", donor_ID_vect, "\n", sep = "\t" )
cat( "0.0. resim_conditionals_drivers.R main: length( donor_ID_vect ) = ", length( donor_ID_vect ), "\n", sep = "\t" )
# readLines("stdin",n=1)

cat( "0.0. resim_conditionals_drivers.R main: donor_age_vect = ", donor_age_vect, "\n", sep = "\t" )
cat( "0.0. resim_conditionals_drivers.R main: length( donor_age_vect ) = ", length( donor_age_vect ), "\n", sep = "\t" )
# readLines("stdin",n=1)

cat( "0.0. resim_conditionals_drivers.R main: donor_n_genomes_vect = ", donor_n_genomes_vect, "\n", sep = "\t" )
cat( "0.0. resim_conditionals_drivers.R main: length( donor_n_genomes_vect ) = ", length( donor_n_genomes_vect ), "\n", sep = "\t" )
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


#Set time points
stat_time=c(0.2,0.4,0.6,0.8)
run <- round(runif( n=1, min= 1e6, max = 2e6), digits = 0)
# my_working_directory = ("/lustre/scratch117/casm/team154/em16/simulations/rsimpop_drivers/rsimpop/")
# my_working_directory = output_path

##Source functions
source(paste0(source_directory,"/functions/filters.parallel_functions.R"))
source(paste0(source_directory,"/functions/tree_functions.R"))
source(paste0(source_directory,"/functions/modelling_functions.R"))
source(paste0(source_directory,"/functions/targeted_analysis_functions.R"))
setwd("/lustre/scratch119/casm/team154pc/ms56/fetal_HSC/treemut"); source("treemut.R");
setwd(output_path)

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
# Order donors w.r.t. age;
# (so that dats for younger donors is generated first)
#

age_order_index_vect = order( donor_age_vect )

donor_ID_vect = donor_ID_vect[ age_order_index_vect ]
donor_age_vect = donor_age_vect[ age_order_index_vect ]
donor_n_genomes_vect = donor_n_genomes_vect[ age_order_index_vect ]

cat( "1.0. resim_conditionals_drivers.R main: donor_ID_vect = ", donor_ID_vect, "\n", sep = "\t" )
cat( "1.0. resim_conditionals_drivers.R main: length( donor_ID_vect ) = ", length( donor_ID_vect ), "\n", sep = "\t" )
# readLines("stdin",n=1)

cat( "1.0. resim_conditionals_drivers.R main: donor_age_vect = ", donor_age_vect, "\n", sep = "\t" )
cat( "1.0. resim_conditionals_drivers.R main: length( donor_age_vect ) = ", length( donor_age_vect ), "\n", sep = "\t" )
# readLines("stdin",n=1)

cat( "1.0. resim_conditionals_drivers.R main: donor_n_genomes_vect = ", donor_n_genomes_vect, "\n", sep = "\t" )
cat( "1.0. resim_conditionals_drivers.R main: length( donor_n_genomes_vect ) = ", length( donor_n_genomes_vect ), "\n", sep = "\t" )
# readLines("stdin",n=1)

#
# Create vectors of (donor-specific) directory paths;
#

donor_output_path_vect = sapply( donor_ID_vect, function(x){ paste( output_path, "/output_", x, "/", sep="" ) } )
resim_output_path_vect = sapply( donor_output_path_vect, function(x){ paste( x, resim_output_directory_name, "/", sep="" ) } )
resim_outfile_vect = sapply( resim_output_path_vect, function(x){ paste( x, "/resim_stats_post_obs_", post_obs_index, "_job_index_", resim_job_index, "_run_", run, ".txt", sep = "" ) } )

####################################################################################################
#
# Create vectors of column names;
# Create data frame "resim_stats_table";
#

# nsim_per_job = 20
nsim_per_job = n_sims_per_job

fitness_threshold = 0.05
number_params=3
cols_params= c("number_drivers_per_year","gamma_shape","gamma_rate")
number_stats=32
cols_stats= c("mean_muts_sim", "n_singletons_1","clade_size_a_1","clade_size_b_1","clade_size_c_1","clade_size_d_1","clade_size_e_1","n_singletons_2","clade_size_a_2","clade_size_b_2","clade_size_c_2","clade_size_d_2","clade_size_e_2", "n_singletons_3", "clade_size_a_3","clade_size_b_3","clade_size_c_3","clade_size_d_3","clade_size_e_3","n_singletons_4", "clade_size_a_4","clade_size_b_4","clade_size_c_4","clade_size_d_4","clade_size_e_4", "ltt_1", "ltt_2", "ltt_3", "ltt_4", "coal_1", "coal_2", "coal_3")

cat( "1.1. resim_conditionals_drivers.R main: length( cols_stats ) = ", length( cols_stats ), "\n", sep = "\t" )
# readLines("stdin",n=1)
	  
n_top_drivers = 5
top_drivers_col_names_vect = sapply( 1:n_top_drivers, function(x){ paste( "driver_fitness_", x, sep="" ) } )

n_top_clades = 5
top_clades_col_names_vect = sapply( 1:n_top_clades, function(x){ paste( "clade_fitness_", x, sep="" ) } )

driver_count_col_names_vect = c( "shared_drivers", "private_drivers", "total_drivers" )

hidden_variable_col_names_vect = c( top_drivers_col_names_vect, top_clades_col_names_vect, driver_count_col_names_vect )

resim_stats_col_names_vect = c( "post_obs_index", cols_stats, hidden_variable_col_names_vect )
resim_stats_table = as.data.frame( matrix( NA, nrow = 1, ncol = length( resim_stats_col_names_vect ) ) )
names( resim_stats_table ) <- resim_stats_col_names_vect

resim_stats_table[[ "post_obs_index" ]][ 1 ] = post_obs_index

####################################################################################################
#
# Read parameter vector from file;
#

posterior_table = read.table( file = posterior_sample_file, sep = "\t", header=TRUE, stringsAsFactors = FALSE )
 
number_drivers_per_year = posterior_table[[ "number_drivers_per_year" ]][ post_obs_index ]
gamma_shape = posterior_table[[ "gamma_shape" ]][ post_obs_index ]
gamma_rate = posterior_table[[ "gamma_rate" ]][ post_obs_index ]

####################################################################################################
#
# Run simulations;
# Fill data frame "resim_stats_table";
#

#target_pop_size = 1e5
#final_division_rate = 0.5/365

for( i in 1:nsim_per_job) 
{
  set.seed(run+i)
  
cat( "2.0. resim_conditionals_drivers.R main: i = ", i, "\n", sep = "\t" )
# readLines("stdin",n=1)
  
  param_ID = paste(number_drivers_per_year, gamma_shape, gamma_rate, sep = "_")
  parameters = list(number_drivers_per_year = number_drivers_per_year, gamma_shape = gamma_shape, gamma_rate = gamma_rate)
  print(parameters)
  
  ##Function to generate gamma distribution based fitness
  genGammaFitness=function(fitness_threshold,shape,rate){
    function() rtrunc(n=1,a=fitness_threshold, b=Inf,"gamma",shape=shape,rate=rate)
  }
  fitnessGammaFn=genGammaFitness(fitness_threshold=fitness_threshold,shape = gamma_shape, rate=gamma_rate)
  
  ####################################################################################################
  
  for( donor_index in 1:length( donor_ID_vect ) )
  {
	  donor_ID = donor_ID_vect[ donor_index ]
	  donor_age = donor_age_vect[ donor_index ]
	  donor_n_genomes = donor_n_genomes_vect[ donor_index ]
	  
	  donor_output_path = donor_output_path_vect[ donor_index ]
	  resim_output_path = resim_output_path_vect[ donor_index ]
	  resim_outfile = resim_outfile_vect[ donor_index ]
	  
  ####################################################################################################
	  
	  if( donor_index == 1 )
	  {
		  ### Run the simulation for donor
		  # dps = run_driver_process_sim( 0.1, 0.5/(365), target_pop_size = 1e5, nyears = donor_age, fitness= fitnessGammaFn, drivers_per_year = number_drivers_per_year )
		  dps = run_driver_process_sim_KEVIN( 0.1, 0.5/(365), target_pop_size = 1e5, nyears = donor_age, fitness= fitnessGammaFn, drivers_per_year = number_drivers_per_year )
		  
	  }else
	  {
		  
		  ### Continue the simulation
		  # dps = continue_driver_process_sim( dps, nyears = donor_age, fitness= fitnessGammaFn )
		  dps = continue_driver_process_sim_KEVIN( dps, nyears = donor_age, fitness= fitnessGammaFn )
		  
	  }
	  
  ####################################################################################################
	  
	  ###Get subsampled tree
	  # dpst_KX001=get_subsampled_tree(dps_KX001,KX001_WGS_sequenced)
	  dpst = get_subsampled_tree( dps, donor_n_genomes )
	  
	  # tree=get_elapsed_time_tree(dpst_KX001,mutrateperdivision = 1,backgroundrate = 16/365)
	  tree = get_elapsed_time_tree( dpst, mutrateperdivision = 1, backgroundrate = 16/365 )
	  
	  mean_muts_sim=mean(get_mut_burden(tree))
	  
	  tree <- make.ultrametric.tree(tree)
	  
  ####################################################################################################
	  
	  #
	  # Get summary statistics;
	  #
	  
	  summary_stats <- get_summary_stats(tree,time_points=stat_time)
	  colnames(summary_stats) <- c("n_singletons_1","clade_size_a_1","clade_size_b_1","clade_size_c_1","clade_size_d_1","clade_size_e_1","n_singletons_2","clade_size_a_2","clade_size_b_2","clade_size_c_2","clade_size_d_2","clade_size_e_2", "n_singletons_3", "clade_size_a_3","clade_size_b_3","clade_size_c_3","clade_size_d_3","clade_size_e_3","n_singletons_4", "clade_size_a_4","clade_size_b_4","clade_size_c_4","clade_size_d_4","clade_size_e_4")
	  ltt_sim <- get_ltt(tree,time_points=stat_time); names(ltt_sim)=paste("ltt",1:length(ltt_sim),sep="_")
	  coal_sim <- get_coalescences(ltt_sim);names(coal_sim)=paste("coal",1:length(coal_sim),sep="_")
	  stats_sim = cbind(mean_muts_sim,summary_stats,t(ltt_sim),t(coal_sim))
	  
	  # sim_table[ i, cols_stats ] <- unlist(stats_sim)
	  resim_stats_table[ 1, cols_stats ] <- unlist(stats_sim)
	  
cat( "2.1. resim_conditionals_drivers.R main: dim( summary_stats ) = ", dim( summary_stats ), "\n", sep = "\t" )
cat( "2.1. resim_conditionals_drivers.R main: dim( stats_sim ) = ", dim( stats_sim ), "\n", sep = "\t" )
cat( "2.1. resim_conditionals_drivers.R main: length( unlist(stats_sim) ) = ", length( unlist(stats_sim) ), "\n", sep = "\t" )
# readLines("stdin",n=1)
	  
	  #
	  # Get hidden variables;
	  #
	  
	  top_drivers_table = get_top_drivers( dps, top.n=n_top_drivers )
	  top_clades_table = get_top_clades( dps, top.n=n_top_clades, min.pop=1 )
	  driver_count_list = get_driver_counts( dps )
	  
	  print( top_drivers_table )
	  
cat( "2.2. resim_conditionals_drivers.R main: dim( top_drivers_table ) = ", dim( top_drivers_table ), "\n", sep = "\t" )
# readLines("stdin",n=1)
	  
	  print( top_clades_table )
	  
cat( "2.2. resim_conditionals_drivers.R main: dim( top_clades_table ) = ", dim( top_clades_table ), "\n", sep = "\t" )
# readLines("stdin",n=1)
	  
	  print( driver_count_list )
	  
cat( "2.2. resim_conditionals_drivers.R main: length( driver_count_list ) = ", length( driver_count_list ), "\n", sep = "\t" )
# readLines("stdin",n=1)
	  
	  resim_stats_table[ 1, top_drivers_col_names_vect ] <- top_drivers_table[[ "fitness" ]]
	  
	  resim_stats_table[ 1, top_clades_col_names_vect ] <- top_clades_table[[ "fitness" ]]
	  
	  resim_stats_table[ 1, driver_count_col_names_vect ] <- unlist( driver_count_list )
	  
  ####################################################################################################
	  
	  # resim_outfile = resim_outfile_vect[ donor_index ]
	  
	  if( i == 1 )
	  {
		  write.table( resim_stats_table[ 1, ], file = resim_outfile , append = FALSE , quote = FALSE , sep = "\t" , row.names = FALSE , col.names = TRUE )
		  # write.table( resim_stats_table[ 1, ], file = resim_outfile , append = TRUE , quote = FALSE , sep = "\t" , row.names = FALSE , col.names = FALSE )
		  
	  }else
	  {
		  # write.table( resim_stats_table[ 1, ], file = resim_outfile , append = FALSE , quote = FALSE , sep = "\t" , row.names = FALSE , col.names = TRUE )
		  write.table( resim_stats_table[ 1, ], file = resim_outfile , append = TRUE , quote = FALSE , sep = "\t" , row.names = FALSE , col.names = FALSE )
		  
	  }
	  
  } # for( donor_index in 1:length( donor_ID_vect ) )
  
} # for( i in 1:nsim_per_job )

####################################################################################################

cat( "done!", "\n" )

