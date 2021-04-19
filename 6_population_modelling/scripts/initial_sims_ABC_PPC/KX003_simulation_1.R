## Run on farm5
args = commandArgs(TRUE)

##SET param_ID AND FILEPATHS
RUNID = as.character(args[1])
output_directory = toString(args[2])
run <- round(runif(1, min= 1e6, max = 2e6), digits = 0)

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
Iteration = "KX003_5"
stat_time=c(0.25,0.50,0.75)
my_working_directory = (paste0("/lustre/scratch117/casm/team154/em16/simulations/rsimpop2/",Iteration))

##Initialise the environment
source("/lustre/scratch117/casm/team154/em16/simulations/rsimpop2/setup.R")
source("/lustre/scratch117/casm/team154/em16/R_scripts/functions/summary_stats_functions.R")
source("/lustre/scratch117/casm/team154/em16/R_scripts/functions/foetal.filters.parallel.R")
source("/lustre/scratch117/casm/team154/em16/R_scripts/functions/phytools_scripts.R")
source("/lustre/scratch117/casm/team154/em16/R_scripts/functions/targeted_analysis_functions.R")
source("/lustre/scratch117/casm/team154/em16/R_scripts/functions/Pop_size_estimation_functions.R")
setwd("/lustre/scratch119/casm/team154pc/ms56/fetal_HSC/treemut"); source("treemut.R");setwd(my_working_directory)

#Specify inputs - need to edit for each individual with cut offs based on phylodyn plots
WGS_sequenced = 328
nsim_per_job = 200
age = 81
midlife = 20:50
latelife = 61:age
change1=15:19
change2=51:60
number_params=3
cols_params= c("target_pop","midlife_pop_change","latelife_pop_change")
number_stats=12
cols_stats= c("mean_muts_sim", "n_singletons_1","largest_clade_1","n_singletons_2", "largest_clade_2", "n_singletons_3", "largest_clade_3", "ltt_1", "ltt_2", "ltt_3", "coal_1", "coal_2")

##Set up output matrices
params.sim=matrix(NA,nrow=nsim_per_job,ncol=number_params)
colnames(params.sim)<-cols_params
stats.sim=matrix(NA,nrow=nsim_per_job,ncol=number_stats)
colnames(stats.sim)<-cols_stats 

##Setup the trajectory
for(i in 1:nsim_per_job) {
  set.seed(run+i)
  target_pop_size= round(runif(1, min= 1e3, max = 2.5e5), digits = 0)
  midlife_pop_change = round(runif(1, min=0.01, max=1.2), digits = 4)
  latelife_pop_change = round(runif(1, min = 0.5, max = 8), digits = 4)
  dif1 = ((target_pop_size*midlife_pop_change)-target_pop_size)/5
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
  
  ##Run the sim
  sp=run_neutral_trajectory(NULL,0.5,trajectory)
  
  ## Visualise the trajectory
  #pdf(paste0("pdfs/",ID, "_", param_ID, "_trajectory.pdf"))
  #plot(sp,xlim=c(0,100))
  #lines(trajectory$ts/365,trajectory$target_pop_size,col="red")
  #legend("topright",c("Target","Actual"),col=c("red","black"),lwd=1)
  #dev.off()
  
  ## Get the tree of extant lineages
  fulltree=get_tree_from_simpop(sp)
  
  ## Subsample
  st=get_subsampled_tree(fulltree,WGS_sequenced)
  
  ## Get elapsed time tree
  tree = get_elapsed_time_tree(st, mutrateperdivision = 1, backgroundrate = 15/365)
  
  mean_muts_sim=mean(get_mut_burden(tree))
  
  tree <- make.ultrametric.tree(tree)
  
  ## The plot true-time version of tree
  #plot_tree(get_elapsed_time_tree(st),cex.label = 0)
  
  #pdf(paste0("pdfs/",ID, "_", param_ID, "_tree.pdf"), width =30)
  #plot_tree(tree,cex.label = 0)
    
  #dev.off()
  
  #Define functions needed for summary stats
  get_summary_stats=function(tree,time_points) {
    nodeheights <- nodeHeights(tree)
    sumstats_each_timepoint=lapply(time_points, function(time_point) {
      clades_post_cutoff=tree$edge[,2][nodeheights[,1] < time_point & !nodeheights[,2] < time_point]
      clade_sizes=sapply(clades_post_cutoff,function(node) {length(getTips(tree,node))})
      n_singletons=sum(clade_sizes==1)
      largest_clade=max(clade_sizes)
      return(data.frame(n_singletons=n_singletons,largest_clade=largest_clade))
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
  colnames(summary_stats) <- c("n_singletons_1","largest_clade_1","n_singletons_2", "largest_clade_2", "n_singletons_3", "largest_clade_3")
  ltt_sim <- get_ltt(tree,time_points=stat_time); names(ltt_sim)=paste("ltt",1:length(ltt_sim),sep="_")
  coal_sim <- get_coalescences(ltt_sim);names(coal_sim)=paste("coals",1:length(coal_sim),sep="_")
  stats_sim = cbind(mean_muts_sim,summary_stats,t(ltt_sim),t(coal_sim))
  
  params.sim[i,] <- unlist(parameters)
  stats.sim[i,]<- unlist(stats_sim)
  
  

  #plot_tree(get_elapsed_time_tree(st),cex.label = 0)
  
    #pdf(paste0("pdfs/",ID, "_", param_ID, "_tree.pdf"), height =30)
    #plot(tree, show.tip.label = FALSE ); abline(v=stat_time,col="red")
    #dev.off()

}

save(stats.sim, file = paste0("summary_stats/summary_stats_run_", run))
save(params.sim, file = paste0("parameters/parameters_run_", run))
