## Run interactively on farm

#Author: Emily Mitchell

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
ID = "KX003"
Iteration = "KX003_5"
Run_ID = "KX003_5_01"
my_working_directory = (paste0("/lustre/scratch117/casm/team154/em16/simulations/rsimpop2/",Iteration))
setwd = my_working_directory
tree_path = paste0("/lustre/scratch117/casm/team154/em16/KX003/KX003_5/trees_vaf_cut/tree_",Run_ID,"_standard_rho01_final.tree")
stat_time=c(0.25,0.50,0.75)

##Initialise the environment
source("/lustre/scratch117/casm/team154/em16/R_scripts/functions/filters_parallel_functions.R")
source("/lustre/scratch117/casm/team154/em16/R_scripts/functions/tree_functions.R")
setwd("/lustre/scratch119/casm/team154pc/ms56/fetal_HSC/treemut"); source("treemut.R");setwd(my_working_directory)

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

#Get the data tree in shape and define the summary statistics 
system("mkdir pdfs")
tree.data=read.tree(tree_path) #Import the observed data tree
mean_muts_data=mean(get_mut_burden(tree.data)) #Calculate the mean mutations per cell stat before make tree ultrametric
tree.data<-make.ultrametric.tree(tree.data) #Make ultrametric
pdf(paste0("pdfs/",ID, "_obs_tree.pdf"), height =30)
plot(tree.data, show.tip.label = FALSE ); abline(v=stat_time,col="red")
dev.off()

summary_stats <- get_summary_stats(tree.data,time_points=stat_time)
colnames(summary_stats) <- c("n_singletons_1","largest_clade_1","n_singletons_2", "largest_clade_2", "n_singletons_3", "largest_clade_3")
ltt_data <- get_ltt(tree.data,time_points=stat_time); names(ltt_data)=paste("ltt",1:length(ltt_data),sep="_")
coal_data <- get_coalescences(ltt_data);names(coal_data)=paste("coals",1:length(coal_data),sep="_")
stats_data = cbind(mean_muts_data,summary_stats,t(ltt_data),t(coal_data))
assign(paste0(ID,".summary_stats.data"),stats_data)

#Define list of simulation tree and parameter files
parameter_run_files=list.files(paste0(my_working_directory,"/parameters/"),pattern="parameters")
stats_run_files=list.files(paste0(my_working_directory,"/summary_stats/"),pattern="summary_stats")

#Create the full parameters matrix
load(paste0(my_working_directory,"/parameters/",parameter_run_files[1]))
parameters.sim <- unlist(params.sim)
for(i in 2:length(parameter_run_files)) {
  load(paste0(my_working_directory,"/parameters/",parameter_run_files[i]))
  parameters.sim  <- rbind(parameters.sim, unlist(params.sim))
}

#Create the full summary_stats matrix
load(paste0(my_working_directory,"/summary_stats/",stats_run_files[1]))
summary_stats.sim <- unlist(stats.sim)
for(i in 2:length(stats_run_files)) {
  load(paste0(my_working_directory,"/summary_stats/",stats_run_files[i]))
  summary_stats.sim  <- rbind(summary_stats.sim, unlist(stats.sim))
}

#Assign the key objects to objects incorporating the model name
assign(paste0(ID,".parameters.sim"),parameters.sim)
assign(paste0(ID,".summary_stats.sim"),summary_stats.sim)

save("KX003.parameters.sim",
     "KX003.summary_stats.sim",
     "KX003.summary_stats.data",
     file = paste0(ID,"_tables"))

#NOW DO THE ABC
#First need to select the summary stats that will be used
all_stats = 1:ncol(summary_stats.sim)
ltt_stats = grep("ltt",colnames(summary_stats.sim))
coals_stats = grep("coals",colnames(summary_stats.sim))
clade_stats = grep("clade",colnames(summary_stats.sim))
singleton_stats = grep("singletons",colnames(summary_stats.sim))

sumstats_to_include = c(ltt_stats)

#Define which model is being used and assign the objects approproiately
model_for_abc = ID
param = get(paste0(ID,".parameters.sim"))
sumstat=get(paste0(ID,".summary_stats.sim"))
sumstat=sumstat[,sumstats_to_include]
target = get(paste0(ID,".summary_stats.data"))
target=target[sumstats_to_include]

tol = 0.01

assign(paste0("rej.",model_for_abc),abc(target = target, param = param, sumstat = sumstat, tol = tol, method = "rejection"))
assign(paste0("param.abc.",model_for_abc),param)

data = as.data.frame(rej.KX003$unadj.values); data$type="data"
prior = as.data.frame(param); prior$type="prior"
data_and_prior = rbind(data,prior)

pdf(paste0("pdfs/",ID, "_data_prior_target_pop.pdf"))
ggplot(data_and_prior, aes(x=target_pop, col=factor(type))) + geom_density() + theme_classic()
dev.off()
pdf(paste0("pdfs/",ID, "_data_prior_latelife.pdf"))
ggplot(data_and_prior, aes(x=latelife_pop_change, col=factor(type))) + geom_density() + theme_classic()
dev.off()
pdf(paste0("pdfs/",ID, "_data_prior_midlife.pdf"))
ggplot(data_and_prior, aes(x=midlife_pop_change, col=factor(type))) + geom_density() + theme_classic()
dev.off()

#Calculate and print quantile results
latelife_result_95CI=quantile(rej.KX003$unadj.values[,"latelife_pop_change"],probs=c(0.025,0.975))
latelife_result_80CI=quantile(rej.KX003$unadj.values[,"latelife_pop_change"],probs=c(0.1,0.9))
latelife_result_95CI
latelife_result_80CI

midlife_result_95CI=quantile(rej.KX003$unadj.values[,"midlife_pop_change"],probs=c(0.025,0.975))
midlife_result_80CI=quantile(rej.KX003$unadj.values[,"midlife_pop_change"],probs=c(0.1,0.9))
midlife_result_95CI
midlife_result_80CI

target_pop_result_95CI=quantile(rej.KX003$unadj.values[,"target_pop"],probs=c(0.025,0.975))
target_pop_result_80CI=quantile(rej.KX003$unadj.values[,"target_pop"],probs=c(0.1,0.9))
target_pop_result_95CI
target_pop_result_80CI

#Save accepted parameters as 5 separate trajectory input files for re-running simulations
traj1 <- rej.KX003$unadj.values[1:200,]
traj2 <- rej.KX003$unadj.values[201:400,]
traj3 <- rej.KX003$unadj.values[401:600,]
traj4 <- rej.KX003$unadj.values[601:800,]
traj5 <- rej.KX003$unadj.values[801:1000,]

save("traj1","traj2","traj3","traj4","traj5",file=paste0(ID,"_traj_for_sim_2"))