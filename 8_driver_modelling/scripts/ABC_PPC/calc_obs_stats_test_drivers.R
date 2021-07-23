
####################################################################################################
#
# "Round 1" simulations are the simulations used to perform ABC sampling.
# "Round 2" simulations ("resimulations") are the simulations used to perform PPC calculations.
#
# As part of the ABC computation at the end of "round 1",
# the observed values of the summary satatistics are computed (for each individual)
# by the script: run_ABC_drivers.R;
# 
# Only when the observed values of the summary satatistics have not been computed for all individuals
# (by running script: resim_conditionals_drivers.R);
# do we need to run the present script (compute_conditionals_drivers.R) to complete this task;
#
# The observed values of the summary satatistics are computed (for each individual)
# are needed as inputs for the script: compute_conditionals_drivers.R;
#
# The present script (compute_conditionals_drivers.R) also collects and pools all these samples from the conditional distributions
# to obtain a stratified sample from the posterior predictive distribution for each summary statistic;
#
####################################################################################################

## Run on farm5
args = commandArgs(TRUE)
args

source_directory = toString(args[1])
tree_file = toString(args[2])
output_directory = toString(args[3])
obs_stats_file_name = toString(args[4])
patient_ID = toString(args[5])

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

cols_params= c("number_drivers_per_year","gamma_shape","gamma_rate")
cols_stats= c("mean_muts_sim", "n_singletons_1","clade_size_a_1","clade_size_b_1","clade_size_c_1","clade_size_d_1","clade_size_e_1","n_singletons_2","clade_size_a_2","clade_size_b_2","clade_size_c_2","clade_size_d_2","clade_size_e_2", "n_singletons_3", "clade_size_a_3","clade_size_b_3","clade_size_c_3","clade_size_d_3","clade_size_e_3","n_singletons_4", "clade_size_a_4","clade_size_b_4","clade_size_c_4","clade_size_d_4","clade_size_e_4", "ltt_1", "ltt_2", "ltt_3", "ltt_4", "coal_1", "coal_2", "coal_3")

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
# assign(paste0(ID,".summary_stats.data"),stats_data)

####################################################################################################
#
# Read summary stats from observed data to file;
#

obs_stats_table = as.data.frame( matrix( NA, nrow = 1, ncol = length( cols_stats ) ) )
names( obs_stats_table ) <- cols_stats

obs_stats_table[ 1, cols_stats ] <- unlist(stats_data)

obs_stats_file = paste0( output_directory, "/", obs_stats_file_name )
write.table( obs_stats_table, file = obs_stats_file , append = FALSE , quote = FALSE , sep = "\t" , row.names = FALSE , col.names = TRUE )

####################################################################################################

cat( "done!", "\n" )

