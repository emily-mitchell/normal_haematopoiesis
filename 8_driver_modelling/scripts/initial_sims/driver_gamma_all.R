###Driver modelling using rsimpop

## Run on farm5
args = commandArgs(TRUE)

##SET param_ID AND FILEPATHS
RUNID = as.character(args[1])
output_directory = toString(args[2])

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

#Set working directory
Iteration = "all"
stat_time=c(0.2,0.4,0.6,0.8)
run <- round(runif(1, min= 1e6, max = 2e6), digits = 0)
my_working_directory = ("/lustre/scratch117/casm/team154/em16/simulations/rsimpop_drivers/rsimpop/")

##Source functions
source("/lustre/scratch117/casm/team154/em16/R_scripts/functions/summary_stats_functions.R")
source("/lustre/scratch117/casm/team154/em16/R_scripts/functions/foetal.filters.parallel.R")
source("/lustre/scratch117/casm/team154/em16/R_scripts/functions/phytools_scripts.R")
source("/lustre/scratch117/casm/team154/em16/R_scripts/functions/targeted_analysis_functions.R")
source("/lustre/scratch117/casm/team154/em16/R_scripts/functions/Pop_size_estimation_functions.R")
source("/lustre/scratch117/casm/team154/em16/R_scripts/functions/tree_plot_fix.R")
setwd("/lustre/scratch119/casm/team154pc/ms56/fetal_HSC/treemut"); source("treemut.R");setwd(my_working_directory)

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

#Specify inputs 
KX001_WGS_sequenced = 408
KX001_age = 29

KX002_WGS_sequenced = 380
KX002_age = 38

SX001_WGS_sequenced = 363
SX001_age = 48

AX001_WGS_sequenced = 368
AX001_age = 63

KX007_WGS_sequenced = 315
KX007_age = 75

KX008_WGS_sequenced = 367
KX008_age = 76

KX004_WGS_sequenced = 451
KX004_age = 77

KX003_WGS_sequenced = 328
KX003_age = 81

nsim_per_job = 20

fitness_threshold = 0.05
number_params=3
cols_params= c("number_drivers_per_year","gamma_shape","gamma_rate")
number_stats=32
cols_stats= c("mean_muts_sim", "n_singletons_1","clade_size_a_1","clade_size_b_1","clade_size_c_1","clade_size_d_1","clade_size_e_1","n_singletons_2","clade_size_a_2","clade_size_b_2","clade_size_c_2","clade_size_d_2","clade_size_e_2", "n_singletons_3", "clade_size_a_3","clade_size_b_3","clade_size_c_3","clade_size_d_3","clade_size_e_3","n_singletons_4", "clade_size_a_4","clade_size_b_4","clade_size_c_4","clade_size_d_4","clade_size_e_4", "ltt_1", "ltt_2", "ltt_3", "ltt_4", "coal_1", "coal_2", "coal_3")

##Set up output matrices
KX001.params.sim=matrix(NA,nrow=nsim_per_job,ncol=number_params)
colnames(KX001.params.sim)<-cols_params
KX001.stats.sim=matrix(NA,nrow=nsim_per_job,ncol=number_stats)
colnames(KX001.stats.sim)<-cols_stats 

KX002.params.sim=matrix(NA,nrow=nsim_per_job,ncol=number_params)
colnames(KX002.params.sim)<-cols_params
KX002.stats.sim=matrix(NA,nrow=nsim_per_job,ncol=number_stats)
colnames(KX002.stats.sim)<-cols_stats 

SX001.params.sim=matrix(NA,nrow=nsim_per_job,ncol=number_params)
colnames(SX001.params.sim)<-cols_params
SX001.stats.sim=matrix(NA,nrow=nsim_per_job,ncol=number_stats)
colnames(SX001.stats.sim)<-cols_stats 

AX001.params.sim=matrix(NA,nrow=nsim_per_job,ncol=number_params)
colnames(AX001.params.sim)<-cols_params
AX001.stats.sim=matrix(NA,nrow=nsim_per_job,ncol=number_stats)
colnames(AX001.stats.sim)<-cols_stats 

KX007.params.sim=matrix(NA,nrow=nsim_per_job,ncol=number_params)
colnames(KX007.params.sim)<-cols_params
KX007.stats.sim=matrix(NA,nrow=nsim_per_job,ncol=number_stats)
colnames(KX007.stats.sim)<-cols_stats 

KX008.params.sim=matrix(NA,nrow=nsim_per_job,ncol=number_params)
colnames(KX008.params.sim)<-cols_params
KX008.stats.sim=matrix(NA,nrow=nsim_per_job,ncol=number_stats)
colnames(KX008.stats.sim)<-cols_stats 

KX004.params.sim=matrix(NA,nrow=nsim_per_job,ncol=number_params)
colnames(KX004.params.sim)<-cols_params
KX004.stats.sim=matrix(NA,nrow=nsim_per_job,ncol=number_stats)
colnames(KX004.stats.sim)<-cols_stats 

KX003.params.sim=matrix(NA,nrow=nsim_per_job,ncol=number_params)
colnames(KX003.params.sim)<-cols_params
KX003.stats.sim=matrix(NA,nrow=nsim_per_job,ncol=number_stats)
colnames(KX003.stats.sim)<-cols_stats 

for(i in 1:nsim_per_job) {
  set.seed(run+i)
  gamma_shape = round(runif(1, min = 0.1, max =2), digits = 2)
  gamma_rate = round(runif(1, min = 5, max = 120), digits = 2)
  number_drivers_per_year = round(runif(1, min =1, max = 200), digits = 0)
  param_ID = paste(number_drivers_per_year, gamma_shape, gamma_rate, sep = "_")
  parameters = list(number_drivers_per_year = number_drivers_per_year, gamma_shape = gamma_shape, gamma_rate = gamma_rate)
  print(parameters)
  
  ##Function to generate gamma distribution based fitness
  genGammaFitness=function(fitness_threshold,shape,rate){
    function() rtrunc(n=1,a=fitness_threshold, b=Inf,"gamma",shape=shape,rate=rate)
  }
  fitnessGammaFn=genGammaFitness(fitness_threshold=fitness_threshold,shape = gamma_shape, rate=gamma_rate)

####################################################################################################################################################
# KX001
  
### Run the simulation for KX001
dps_KX001=run_driver_process_sim(0.1, final_division_rate = 0.5/(365),target_pop_size = 1e5, nyears = KX001_age, fitness= fitnessGammaFn, drivers_per_year = number_drivers_per_year)

### Save final per driver counts for KX001
driver_counts <- dps_KX001$cfg$info %>% filter(population>0)
#write.table(driver_counts, paste0("/lustre/scratch117/casm/team154/em16/simulations/rsimpop_drivers/rsimpop/KX001_4/drivers/KX001_4_",param_ID,"_drivers.txt"))

###Get subsampled tree
dpst_KX001=get_subsampled_tree(dps_KX001,KX001_WGS_sequenced)

tree=get_elapsed_time_tree(dpst_KX001,mutrateperdivision = 1,backgroundrate = 16/365)

mean_muts_sim=mean(get_mut_burden(tree))

#pdf(paste0("/lustre/scratch117/casm/team154/em16/simulations/rsimpop_drivers/rsimpop/KX001_4/pdfs/KX001_4_", param_ID, "_tree.pdf"), width =30)
#plot_tree_events2(tree,fmode =1, cex.label = 0)

#dev.off()

tree <- make.ultrametric.tree(tree)

summary_stats <- get_summary_stats(tree,time_points=stat_time)
colnames(summary_stats) <- c("n_singletons_1","clade_size_a_1","clade_size_b_1","clade_size_c_1","clade_size_d_1","clade_size_e_1","n_singletons_2","clade_size_a_2","clade_size_b_2","clade_size_c_2","clade_size_d_2","clade_size_e_2", "n_singletons_3", "clade_size_a_3","clade_size_b_3","clade_size_c_3","clade_size_d_3","clade_size_e_3","n_singletons_4", "clade_size_a_4","clade_size_b_4","clade_size_c_4","clade_size_d_4","clade_size_e_4")
ltt_sim <- get_ltt(tree,time_points=stat_time); names(ltt_sim)=paste("ltt",1:length(ltt_sim),sep="_")
coal_sim <- get_coalescences(ltt_sim);names(coal_sim)=paste("coal",1:length(coal_sim),sep="_")
stats_sim = cbind(mean_muts_sim,summary_stats,t(ltt_sim),t(coal_sim))

KX001.params.sim[i,] <- unlist(parameters)
KX001.stats.sim[i,]<- unlist(stats_sim)




####################################################################################################################################################
# KX002

### Continue the simulation for KX002
dps_KX002 =continue_driver_process_sim(dps_KX001 ,nyears = KX002_age ,fitnessGen = fitnessGammaFn)

### Save final per driver counts for KX002
driver_counts <- dps_KX002$cfg$info %>% filter(population>0)
#write.table(driver_counts, paste0("/lustre/scratch117/casm/team154/em16/simulations/rsimpop_drivers/rsimpop/KX002_2/drivers/KX002_2_",param_ID,"_drivers.txt"))

###Get subsampled tree
dpst_KX002=get_subsampled_tree(dps_KX002,KX002_WGS_sequenced)

tree=get_elapsed_time_tree(dpst_KX002,mutrateperdivision = 1,backgroundrate = 16/365)

mean_muts_sim=mean(get_mut_burden(tree))

#pdf(paste0("/lustre/scratch117/casm/team154/em16/simulations/rsimpop_drivers/rsimpop/KX002_2/pdfs/KX002_2_", param_ID, "_tree.pdf"), width =30)
#plot_tree_events2(tree,fmode =1, cex.label = 0)

#dev.off()

tree <- make.ultrametric.tree(tree)

summary_stats <- get_summary_stats(tree,time_points=stat_time)
colnames(summary_stats) <- c("n_singletons_1","clade_size_a_1","clade_size_b_1","clade_size_c_1","clade_size_d_1","clade_size_e_1","n_singletons_2","clade_size_a_2","clade_size_b_2","clade_size_c_2","clade_size_d_2","clade_size_e_2", "n_singletons_3", "clade_size_a_3","clade_size_b_3","clade_size_c_3","clade_size_d_3","clade_size_e_3","n_singletons_4", "clade_size_a_4","clade_size_b_4","clade_size_c_4","clade_size_d_4","clade_size_e_4")
ltt_sim <- get_ltt(tree,time_points=stat_time); names(ltt_sim)=paste("ltt",1:length(ltt_sim),sep="_")
coal_sim <- get_coalescences(ltt_sim);names(coal_sim)=paste("coal",1:length(coal_sim),sep="_")
stats_sim = cbind(mean_muts_sim,summary_stats,t(ltt_sim),t(coal_sim))

KX002.params.sim[i,] <- unlist(parameters)
KX002.stats.sim[i,]<- unlist(stats_sim)



####################################################################################################################################################
# SX001

### Continue the simulation for SX001
dps_SX001 =continue_driver_process_sim(dps_KX002 ,nyears = SX001_age ,fitnessGen = fitnessGammaFn)

### Save final per driver counts for SX001
driver_counts <- dps_KX002$cfg$info %>% filter(population>0)
#write.table(driver_counts, paste0("/lustre/scratch117/casm/team154/em16/simulations/rsimpop_drivers/rsimpop/SX001_5/drivers/SX001_5_",param_ID,"_drivers.txt"))

###Get subsampled tree
dpst_SX001=get_subsampled_tree(dps_SX001,SX001_WGS_sequenced)

tree=get_elapsed_time_tree(dpst_SX001,mutrateperdivision = 1,backgroundrate = 16/365)

mean_muts_sim=mean(get_mut_burden(tree))

#pdf(paste0("/lustre/scratch117/casm/team154/em16/simulations/rsimpop_drivers/rsimpop/SX001_5/pdfs/SX001_5_", param_ID, "_tree.pdf"), width =30)
#plot_tree_events2(tree,fmode =1, cex.label = 0)

#dev.off()

tree <- make.ultrametric.tree(tree)

summary_stats <- get_summary_stats(tree,time_points=stat_time)
colnames(summary_stats) <- c("n_singletons_1","clade_size_a_1","clade_size_b_1","clade_size_c_1","clade_size_d_1","clade_size_e_1","n_singletons_2","clade_size_a_2","clade_size_b_2","clade_size_c_2","clade_size_d_2","clade_size_e_2", "n_singletons_3", "clade_size_a_3","clade_size_b_3","clade_size_c_3","clade_size_d_3","clade_size_e_3","n_singletons_4", "clade_size_a_4","clade_size_b_4","clade_size_c_4","clade_size_d_4","clade_size_e_4")
ltt_sim <- get_ltt(tree,time_points=stat_time); names(ltt_sim)=paste("ltt",1:length(ltt_sim),sep="_")
coal_sim <- get_coalescences(ltt_sim);names(coal_sim)=paste("coal",1:length(coal_sim),sep="_")
stats_sim = cbind(mean_muts_sim,summary_stats,t(ltt_sim),t(coal_sim))

SX001.params.sim[i,] <- unlist(parameters)
SX001.stats.sim[i,]<- unlist(stats_sim)




####################################################################################################################################################
# AX001

### Continue the simulation for AX001
dps_AX001 =continue_driver_process_sim(dps_SX001 ,nyears = AX001_age ,fitnessGen = fitnessGammaFn)

### Save final per driver counts for AX001
driver_counts <- dps_AX001$cfg$info %>% filter(population>0)
#write.table(driver_counts, paste0("/lustre/scratch117/casm/team154/em16/simulations/rsimpop_drivers/rsimpop/AX001_4/drivers/AX001_4_",param_ID,"_drivers.txt"))

###Get subsampled tree
dpst_AX001=get_subsampled_tree(dps_AX001,AX001_WGS_sequenced)

tree=get_elapsed_time_tree(dpst_AX001,mutrateperdivision = 1,backgroundrate = 16/365)

mean_muts_sim=mean(get_mut_burden(tree))

#pdf(paste0("/lustre/scratch117/casm/team154/em16/simulations/rsimpop_drivers/rsimpop/AX001_4/pdfs/AX001_4_", param_ID, "_tree.pdf"), width =30)
#plot_tree_events2(tree,fmode =1, cex.label = 0)

#dev.off()

tree <- make.ultrametric.tree(tree)

summary_stats <- get_summary_stats(tree,time_points=stat_time)
colnames(summary_stats) <- c("n_singletons_1","clade_size_a_1","clade_size_b_1","clade_size_c_1","clade_size_d_1","clade_size_e_1","n_singletons_2","clade_size_a_2","clade_size_b_2","clade_size_c_2","clade_size_d_2","clade_size_e_2", "n_singletons_3", "clade_size_a_3","clade_size_b_3","clade_size_c_3","clade_size_d_3","clade_size_e_3","n_singletons_4", "clade_size_a_4","clade_size_b_4","clade_size_c_4","clade_size_d_4","clade_size_e_4")
ltt_sim <- get_ltt(tree,time_points=stat_time); names(ltt_sim)=paste("ltt",1:length(ltt_sim),sep="_")
coal_sim <- get_coalescences(ltt_sim);names(coal_sim)=paste("coal",1:length(coal_sim),sep="_")
stats_sim = cbind(mean_muts_sim,summary_stats,t(ltt_sim),t(coal_sim))

AX001.params.sim[i,] <- unlist(parameters)
AX001.stats.sim[i,]<- unlist(stats_sim)



####################################################################################################################################################
# KX007

### Continue the simulation for KX007
dps_KX007 =continue_driver_process_sim(dps_AX001 ,nyears = KX007_age ,fitnessGen = fitnessGammaFn)

### Save final per driver counts for KX007
driver_counts <- dps_KX007$cfg$info %>% filter(population>0)
#write.table(driver_counts, paste0("/lustre/scratch117/casm/team154/em16/simulations/rsimpop_drivers/rsimpop/KX007_2/drivers/KX007_2_",param_ID,"_drivers.txt"))

###Get subsampled tree
dpst_KX007=get_subsampled_tree(dps_KX007,KX007_WGS_sequenced)

tree=get_elapsed_time_tree(dpst_KX007,mutrateperdivision = 1,backgroundrate = 16/365)

mean_muts_sim=mean(get_mut_burden(tree))

#pdf(paste0("/lustre/scratch117/casm/team154/em16/simulations/rsimpop_drivers/rsimpop/KX007_2/pdfs/KX007_2_", param_ID, "_tree.pdf"), width =30)
#plot_tree_events2(tree,fmode =1, cex.label = 0)

#dev.off()

tree <- make.ultrametric.tree(tree)

summary_stats <- get_summary_stats(tree,time_points=stat_time)
colnames(summary_stats) <- c("n_singletons_1","clade_size_a_1","clade_size_b_1","clade_size_c_1","clade_size_d_1","clade_size_e_1","n_singletons_2","clade_size_a_2","clade_size_b_2","clade_size_c_2","clade_size_d_2","clade_size_e_2", "n_singletons_3", "clade_size_a_3","clade_size_b_3","clade_size_c_3","clade_size_d_3","clade_size_e_3","n_singletons_4", "clade_size_a_4","clade_size_b_4","clade_size_c_4","clade_size_d_4","clade_size_e_4")
ltt_sim <- get_ltt(tree,time_points=stat_time); names(ltt_sim)=paste("ltt",1:length(ltt_sim),sep="_")
coal_sim <- get_coalescences(ltt_sim);names(coal_sim)=paste("coal",1:length(coal_sim),sep="_")
stats_sim = cbind(mean_muts_sim,summary_stats,t(ltt_sim),t(coal_sim))

KX007.params.sim[i,] <- unlist(parameters)
KX007.stats.sim[i,]<- unlist(stats_sim)




####################################################################################################################################################
# KX008

### Continue the simulation for KX008
dps_KX008 =continue_driver_process_sim(dps_KX007 ,nyears = KX008_age ,fitnessGen = fitnessGammaFn)

### Save final per driver counts for KX002
driver_counts <- dps_KX008$cfg$info %>% filter(population>0)
#write.table(driver_counts, paste0("/lustre/scratch117/casm/team154/em16/simulations/rsimpop_drivers/rsimpop/KX008_2/drivers/KX008_2_",param_ID,"_drivers.txt"))

###Get subsampled tree
dpst_KX008=get_subsampled_tree(dps_KX008,KX008_WGS_sequenced)

tree=get_elapsed_time_tree(dpst_KX008,mutrateperdivision = 1,backgroundrate = 16/365)

mean_muts_sim=mean(get_mut_burden(tree))

#pdf(paste0("/lustre/scratch117/casm/team154/em16/simulations/rsimpop_drivers/rsimpop/KX008_2/pdfs/KX008_2_", param_ID, "_tree.pdf"), width =30)
#plot_tree_events2(tree,fmode =1, cex.label = 0)

#dev.off()

tree <- make.ultrametric.tree(tree)

summary_stats <- get_summary_stats(tree,time_points=stat_time)
colnames(summary_stats) <- c("n_singletons_1","clade_size_a_1","clade_size_b_1","clade_size_c_1","clade_size_d_1","clade_size_e_1","n_singletons_2","clade_size_a_2","clade_size_b_2","clade_size_c_2","clade_size_d_2","clade_size_e_2", "n_singletons_3", "clade_size_a_3","clade_size_b_3","clade_size_c_3","clade_size_d_3","clade_size_e_3","n_singletons_4", "clade_size_a_4","clade_size_b_4","clade_size_c_4","clade_size_d_4","clade_size_e_4")
ltt_sim <- get_ltt(tree,time_points=stat_time); names(ltt_sim)=paste("ltt",1:length(ltt_sim),sep="_")
coal_sim <- get_coalescences(ltt_sim);names(coal_sim)=paste("coal",1:length(coal_sim),sep="_")
stats_sim = cbind(mean_muts_sim,summary_stats,t(ltt_sim),t(coal_sim))

KX008.params.sim[i,] <- unlist(parameters)
KX008.stats.sim[i,]<- unlist(stats_sim)




####################################################################################################################################################
# KX004

### Continue the simulation for KX004
dps_KX004 =continue_driver_process_sim(dps_KX008 ,nyears = KX004_age ,fitnessGen = fitnessGammaFn)

### Save final per driver counts for KX004
driver_counts <- dps_KX004$cfg$info %>% filter(population>0)
#write.table(driver_counts, paste0("/lustre/scratch117/casm/team154/em16/simulations/rsimpop_drivers/rsimpop/KX004_4/drivers/KX004_4_",param_ID,"_drivers.txt"))

###Get subsampled tree
dpst_KX004=get_subsampled_tree(dps_KX004,KX004_WGS_sequenced)

tree=get_elapsed_time_tree(dpst_KX004,mutrateperdivision = 1,backgroundrate = 16/365)

mean_muts_sim=mean(get_mut_burden(tree))

#pdf(paste0("/lustre/scratch117/casm/team154/em16/simulations/rsimpop_drivers/rsimpop/KX004_4/pdfs/KX004_4_", param_ID, "_tree.pdf"), width =30)
#plot_tree_events2(tree,fmode =1, cex.label = 0)

#dev.off()

tree <- make.ultrametric.tree(tree)

summary_stats <- get_summary_stats(tree,time_points=stat_time)
colnames(summary_stats) <- c("n_singletons_1","clade_size_a_1","clade_size_b_1","clade_size_c_1","clade_size_d_1","clade_size_e_1","n_singletons_2","clade_size_a_2","clade_size_b_2","clade_size_c_2","clade_size_d_2","clade_size_e_2", "n_singletons_3", "clade_size_a_3","clade_size_b_3","clade_size_c_3","clade_size_d_3","clade_size_e_3","n_singletons_4", "clade_size_a_4","clade_size_b_4","clade_size_c_4","clade_size_d_4","clade_size_e_4")
ltt_sim <- get_ltt(tree,time_points=stat_time); names(ltt_sim)=paste("ltt",1:length(ltt_sim),sep="_")
coal_sim <- get_coalescences(ltt_sim);names(coal_sim)=paste("coal",1:length(coal_sim),sep="_")
stats_sim = cbind(mean_muts_sim,summary_stats,t(ltt_sim),t(coal_sim))

KX004.params.sim[i,] <- unlist(parameters)
KX004.stats.sim[i,]<- unlist(stats_sim)




####################################################################################################################################################
# KX003

### Continue the simulation for KX003
dps_KX003 =continue_driver_process_sim(dps_KX004 ,nyears = KX003_age ,fitnessGen = fitnessGammaFn)

### Save final per driver counts for KX003
driver_counts <- dps_KX003$cfg$info %>% filter(population>0)
#write.table(driver_counts, paste0("/lustre/scratch117/casm/team154/em16/simulations/rsimpop_drivers/rsimpop/KX004_4/drivers/KX003_5_",param_ID,"_drivers.txt"))

###Get subsampled tree
dpst_KX003=get_subsampled_tree(dps_KX003,KX003_WGS_sequenced)

tree=get_elapsed_time_tree(dpst_KX003,mutrateperdivision = 1,backgroundrate = 16/365)

mean_muts_sim=mean(get_mut_burden(tree))

#pdf(paste0("/lustre/scratch117/casm/team154/em16/simulations/rsimpop_drivers/rsimpop/KX004_4/pdfs/KX004_4_", param_ID, "_tree.pdf"), width =30)
#plot_tree_events2(tree,fmode =1, cex.label = 0)

#dev.off()

tree <- make.ultrametric.tree(tree)

summary_stats <- get_summary_stats(tree,time_points=stat_time)
colnames(summary_stats) <- c("n_singletons_1","clade_size_a_1","clade_size_b_1","clade_size_c_1","clade_size_d_1","clade_size_e_1","n_singletons_2","clade_size_a_2","clade_size_b_2","clade_size_c_2","clade_size_d_2","clade_size_e_2", "n_singletons_3", "clade_size_a_3","clade_size_b_3","clade_size_c_3","clade_size_d_3","clade_size_e_3","n_singletons_4", "clade_size_a_4","clade_size_b_4","clade_size_c_4","clade_size_d_4","clade_size_e_4")
ltt_sim <- get_ltt(tree,time_points=stat_time); names(ltt_sim)=paste("ltt",1:length(ltt_sim),sep="_")
coal_sim <- get_coalescences(ltt_sim);names(coal_sim)=paste("coal",1:length(coal_sim),sep="_")
stats_sim = cbind(mean_muts_sim,summary_stats,t(ltt_sim),t(coal_sim))

KX003.params.sim[i,] <- unlist(parameters)
KX003.stats.sim[i,]<- unlist(stats_sim)




####################################################################################################################################################


}

save(KX001.stats.sim, file = paste0("/lustre/scratch117/casm/team154/em16/simulations/rsimpop_drivers/rsimpop/KX001_4/summary_stats/summary_stats_run_", run))
save(KX001.params.sim, file = paste0("/lustre/scratch117/casm/team154/em16/simulations/rsimpop_drivers/rsimpop/KX001_4/parameters/parameters_run_", run))

save(KX002.stats.sim, file = paste0("/lustre/scratch117/casm/team154/em16/simulations/rsimpop_drivers/rsimpop/KX002_2/summary_stats/summary_stats_run_", run))
save(KX002.params.sim, file = paste0("/lustre/scratch117/casm/team154/em16/simulations/rsimpop_drivers/rsimpop/KX002_2/parameters/parameters_run_", run))

save(SX001.stats.sim, file = paste0("/lustre/scratch117/casm/team154/em16/simulations/rsimpop_drivers/rsimpop/SX001_5/summary_stats/summary_stats_run_", run))
save(SX001.params.sim, file = paste0("/lustre/scratch117/casm/team154/em16/simulations/rsimpop_drivers/rsimpop/SX001_5/parameters/parameters_run_", run))

save(AX001.stats.sim, file = paste0("/lustre/scratch117/casm/team154/em16/simulations/rsimpop_drivers/rsimpop/AX001_4/summary_stats/summary_stats_run_", run))
save(AX001.params.sim, file = paste0("/lustre/scratch117/casm/team154/em16/simulations/rsimpop_drivers/rsimpop/AX001_4/parameters/parameters_run_", run))

save(KX007.stats.sim, file = paste0("/lustre/scratch117/casm/team154/em16/simulations/rsimpop_drivers/rsimpop/KX007_2/summary_stats/summary_stats_run_", run))
save(KX007.params.sim, file = paste0("/lustre/scratch117/casm/team154/em16/simulations/rsimpop_drivers/rsimpop/KX007_2/parameters/parameters_run_", run))

save(KX008.stats.sim, file = paste0("/lustre/scratch117/casm/team154/em16/simulations/rsimpop_drivers/rsimpop/KX008_2/summary_stats/summary_stats_run_", run))
save(KX008.params.sim, file = paste0("/lustre/scratch117/casm/team154/em16/simulations/rsimpop_drivers/rsimpop/KX008_2/parameters/parameters_run_", run))

save(KX004.stats.sim, file = paste0("/lustre/scratch117/casm/team154/em16/simulations/rsimpop_drivers/rsimpop/KX004_4/summary_stats/summary_stats_run_", run))
save(KX004.params.sim, file = paste0("/lustre/scratch117/casm/team154/em16/simulations/rsimpop_drivers/rsimpop/KX004_4/parameters/parameters_run_", run))

save(KX003.stats.sim, file = paste0("/lustre/scratch117/casm/team154/em16/simulations/rsimpop_drivers/rsimpop/KX003_5/summary_stats/summary_stats_run_", run))
save(KX003.params.sim, file = paste0("/lustre/scratch117/casm/team154/em16/simulations/rsimpop_drivers/rsimpop/KX003_5/parameters/parameters_run_", run))

