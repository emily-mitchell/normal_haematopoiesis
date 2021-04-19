## Run on farm5
args = commandArgs(TRUE)

##SET param_ID AND FILEPATHS
# RUNID = as.character(args[1])
input_file_path = toString(args[1])
output_directory = toString(args[2])
pp_sim_index=as.integer(args[3])
resim_run_index=as.integer(args[4])
n_sims_per_job=as.integer(args[5])
source_directory = toString(args[6])

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

#Set working directory
ID = "KX003"
Iteration = "KX003_5"
stat_time=c(0.2,0.4,0.6,0.8)
my_working_directory = output_directory

##Initialise the environment
setwd("/lustre/scratch117/casm/team154/em16/simulations/rsimpop2/")
source("/lustre/scratch117/casm/team154/em16/simulations/rsimpop2/setup.R")
source("/lustre/scratch117/casm/team154/em16/R_scripts/functions/summary_stats_functions.R")
source("/lustre/scratch117/casm/team154/em16/R_scripts/functions/foetal.filters.parallel.R")
source("/lustre/scratch117/casm/team154/em16/R_scripts/functions/phytools_scripts.R")
source("/lustre/scratch117/casm/team154/em16/R_scripts/functions/targeted_analysis_functions.R")
source("/lustre/scratch117/casm/team154/em16/R_scripts/functions/Pop_size_estimation_functions.R")
setwd("/lustre/scratch119/casm/team154pc/ms56/fetal_HSC/treemut"); source("treemut.R");setwd(my_working_directory)

#Specify inputs
run = resim_run_index
WGS_sequenced = 328
nsim_per_job = n_sims_per_job

age = 81
midlife = 20:50
latelife = 61:age
change1=15:19
change2=51:60

number_params=3
cols_params= c("target_pop","midlife_pop_change","latelife_pop_change")
number_stats=32
cols_stats= c("mean_muts_sim", "n_singletons_1","clade_size_a_1","clade_size_b_1","clade_size_c_1","clade_size_d_1","clade_size_e_1","n_singletons_2","clade_size_a_2","clade_size_b_2","clade_size_c_2","clade_size_d_2","clade_size_e_2", "n_singletons_3", "clade_size_a_3","clade_size_b_3","clade_size_c_3","clade_size_d_3","clade_size_e_3","n_singletons_4", "clade_size_a_4","clade_size_b_4","clade_size_c_4","clade_size_d_4","clade_size_e_4", "ltt_1", "ltt_2", "ltt_3", "ltt_4", "coal_1", "coal_2", "coal_3")

##Set up output matrices
stats.resim=matrix(NA,nrow=nsim_per_job,ncol=number_stats)
colnames(stats.resim)<-cols_stats


##Read in sample from posterior distribution (parameter values);
##Read in sample from posterior predictive distribution (summary stat values);
##Read in observed data set (summary stat values);

load(input_file_path)

pp_param_table <- as.data.frame(KX003.parameters.ppc)

target_pop_size = pp_param_table$target_pop[pp_sim_index]
midlife_pop_change = pp_param_table$midlife_pop_change[pp_sim_index]
latelife_pop_change = pp_param_table$latelife_pop_change[pp_sim_index]

dif1 = ((target_pop_size*midlife_pop_change)-target_pop_size)/15
dif2 = ((target_pop_size*latelife_pop_change)-(target_pop_size*midlife_pop_change))/10
param_ID = paste(target_pop_size, midlife_pop_change, latelife_pop_change, sep = "_")

parameters = list(target_pop= target_pop_size, midlife_pop_change = midlife_pop_change, latelife_pop_change= latelife_pop_change)
print(parameters)

trajectory=data.frame(ts=365*(1:age),target_pop_size=target_pop_size+1*(1:age),division_rate=1/(365))
trajectory$target_pop_size[change1]=round(trajectory$target_pop_size[14]+(dif1*(1:5)),digits=0)
trajectory$target_pop_size[midlife]=round(midlife_pop_change*trajectory$target_pop_size[midlife], digits=0)
trajectory$target_pop_size[change2]=round(trajectory$target_pop_size[50]+(dif2*(1:10)),digits =0)
trajectory$target_pop_size[latelife]=round(latelife_pop_change*trajectory$target_pop_size[latelife], digits=0)
print(head(trajectory))

##Setup the trajectory
for(i in 1:nsim_per_job) {
  set.seed(run+i)
  
  ##Run the sim
  sp=run_neutral_trajectory(NULL,0.5,trajectory)
  
  ## Get the tree of extant lineages
  fulltree=get_tree_from_simpop(sp)
  
  ## Subsample
  st=get_subsampled_tree(fulltree,WGS_sequenced)
  
  ## Get elapsed time tree
  tree = get_elapsed_time_tree(st, mutrateperdivision = 1, backgroundrate = 15/365)
  
  mean_muts_sim=mean(get_mut_burden(tree))
  
  tree <- make.ultrametric.tree(tree)
  
  #Define functions needed for summary stats
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
  
  summary_stats <- get_summary_stats(tree,time_points=stat_time)
  colnames(summary_stats) <- c("n_singletons_1","clade_size_a_1","clade_size_b_1","clade_size_c_1","clade_size_d_1","clade_size_e_1","n_singletons_2","clade_size_a_2","clade_size_b_2","clade_size_c_2","clade_size_d_2","clade_size_e_2", "n_singletons_3", "clade_size_a_3","clade_size_b_3","clade_size_c_3","clade_size_d_3","clade_size_e_3","n_singletons_4", "clade_size_a_4","clade_size_b_4","clade_size_c_4","clade_size_d_4","clade_size_e_4")
  ltt_sim <- get_ltt(tree,time_points=stat_time); names(ltt_sim)=paste("ltt",1:length(ltt_sim),sep="_")
  coal_sim <- get_coalescences(ltt_sim);names(coal_sim)=paste("coal",1:length(coal_sim),sep="_")
  stats_sim = cbind(mean_muts_sim,summary_stats,t(ltt_sim),t(coal_sim))
  
 
  stats.resim[i,]<- unlist(stats_sim)
  
 

} # for(i in 1:nsim_per_job)

save(stats.resim, file = paste0(output_directory,"/summary_stats_pp_sim_", pp_sim_index,"_resim_run_", resim_run_index))
