##BENCH-MARKING OF PHYLOGENY STRUCTURE
##In addition to the bootstrapping which is in a separate script, this script
#1.  Compares the tree to that built by other tree-building algorithms
#2.  Reviews internal consistency of the raw tree genotypes by two methods: (1) testing against perfect assumptions of phylogeny mutations, (2) comparing expected genotypes from the consensus phylogeny vs the actual genotype in the input matrix

#Author: Mike Spencer Chapman

library(stringr)
library(ape)
library(seqinr)
library(ggtree)
library(tidyr)
library(dplyr)
library(ggplot2)
library(plotrix)
library(phangorn)
library(pheatmap)
library(RColorBrewer)

#Define function required later in script for comparing phylogenies of different algorithms
comparePhylo_and_plot=function(tree1,tree2,names){
  plot_comp_tree=function(tree,comp,title,col="red",lwd=1,tree_pos){
    tree_name=deparse(substitute(tree))
    shared_clades=comp[,tree_pos]
    edge_width=sapply(tree$edge[,2],function(node) ifelse(node%in%c(1:length(tree$tip.label),shared_clades),lwd,2*lwd))
    edge_col=sapply(tree$edge[,2],function(node) ifelse(node%in%c(1:length(tree$tip.label),shared_clades),"black",col))
    plot(tree,show.tip.label=F,direction="downwards",edge.color=edge_col,edge.width=edge_width,main=title)
  }
  comp<-compare_nodes(tree1,tree2)
  par(mfrow=c(1,2))
  plot_comp_tree(tree1,comp=comp,title=names[1],tree_pos = 1)
  plot_comp_tree(tree2,comp=comp,title=names[2],tree_pos = 2)
}

compare_nodes=function(x, y) 
{
  tree1 <- deparse(substitute(x))
  tree2 <- deparse(substitute(y))
  n1 <- Ntip(x)
  n2 <- Ntip(y)
  
  key1 <- makeNodeLabel(x, "md5sum")$node.label
  key2 <- makeNodeLabel(y, "md5sum")$node.label
  mk12 <- match(key1, key2)
  mk21 <- match(key2, key1)
  if (any(tmp <- is.na(mk12))) {
    nk <- sum(tmp)
  }
  if (any(tmp <- is.na(mk21))) {
    nk <- sum(tmp)
  }
  nodes1 <- which(!is.na(mk12))
  nodes2 <- mk12[!is.na(mk12)]
  
  NODES <- data.frame(nodes1 + n1,nodes2 + n2)
  names(NODES) <- c(tree1, tree2)
  return(NODES)
}

get_RF_dist=function(tree1,tree2){
  comp<-compare_nodes(tree1,tree2)
  RF_dist=1-(length(comp[,1])/length(unique(tree1$edge[,1])))
  return(RF_dist)
}

#Set working directly & load-up all 'local drive' functions
my_working_directory="/lustre/scratch119/casm/team154pc/ms56/Emily_benchmarking"
setwd(my_working_directory)
R_function_files = list.files("/lustre/scratch119/casm/team154pc/ms56/my_functions",pattern=".R",full.names=TRUE)
treemut_dir="/lustre/scratch119/casm/team154pc/ms56/fetal_HSC/treemut"
sapply(R_function_files[-2],source)
setwd(treemut_dir); source("treemut.R"); setwd(my_working_directory)


all_results=list()

samples=readLines("/lustre/scratch119/casm/team154pc/ms56/lesion_segregation/input_data/EM_samples.txt")

# #Loop to copy over the data from Emily's lustre space
# for(sample in samples) {
#   sample_elements=str_split(sample,pattern="_",simplify=T)
#   annot_muts_file=paste0("/lustre/scratch117/casm/team154/em16/",sample_elements[1],"/",paste0(sample_elements[1:2],collapse="_"),"/trees_vaf_cut/annotated_mut_set_",sample,"_standard_rho01")
#   tree_file=paste0("/lustre/scratch117/casm/team154/em16/",sample_elements[1],"/",paste0(sample_elements[1:2],collapse="_"),"/trees_vaf_cut/tree_",sample,"_standard_rho01.tree")
#   
#   destination_folder=paste0("/lustre/scratch119/casm/team154pc/ms56/Emily_benchmarking/",sample_elements[1],"/vaf_filtered/")
#   copy_command=paste0("cp ",annot_muts_file," ",destination_folder)
#   copy_command2=paste0("cp ",tree_file," ",destination_folder)
#   
#   system(paste0("mkdir -p ",destination_folder))
#   system(copy_command)
#   system(copy_command2)
# }

for(i in 1:length(samples)) {
  this_sample=samples[i]
  #Set file paths
  exp_ID=str_split(this_sample,pattern="_",simplify=T)[1]
  filtering="vaf_filtered"
  run_ID=paste0(str_split(this_sample,pattern="_",simplify=T)[2:3],collapse = "_")
  filtering_ID="standard_rho01"
  my_working_directory=paste0("/lustre/scratch119/casm/team154pc/ms56/Emily_benchmarking/",exp_ID,"/",filtering)
  output_dir=my_working_directory
  
  filtered_muts_file = paste0("annotated_mut_set_",exp_ID,"_",run_ID,"_",filtering_ID)
  tree_file_path=paste0("tree_",exp_ID,"_",run_ID,"_",filtering_ID,".tree")
  
  
  #Load the filtered_muts_file & create the mutation matrix in the format required for scite
  setwd(my_working_directory)
  load(filtered_muts_file)
  gt<-filtered_muts$Genotype_shared_bin
  NV<-as.matrix(filtered_muts$COMB_mats.tree.build$NV);NR<-as.matrix(filtered_muts$COMB_mats.tree.build$NR)
  details<-filtered_muts$COMB_mats.tree.build$mat
  tree<-read.tree(tree_file_path)
  
  
  
  #Creates the "perfect" mutation matrix that would expect from the tree
  gt<-gt[rownames(gt)%in%details$mut_ref,]
  phylogeny_mut_mat=Reduce(rbind,lapply(rownames(gt),function(mut) colnames(gt)%in% getTips(tree,details$node[details$mut_ref==mut])))
  dimnames(phylogeny_mut_mat)=dimnames(gt)
  
  #(b) Comparing pairs of loci for perfect phylogeny assumptions (does not rely on the phylogeny)
  disagreement_score=function(binary_mut_mat) {
    max_mut=5000
    if(dim(binary_mut_mat)[1]>max_mut) {
      binary_mut_mat<-binary_mut_mat[sample(1:nrow(binary_mut_mat),size=max_mut,replace = F),]
    }
    c=matrix(as.logical(binary_mut_mat),ncol=ncol(binary_mut_mat)) #c is a binary matrix with 0.5's rounded up to 1's
    d<-binary_mut_mat;d[d==0.5]<-0;d<-matrix(as.logical(d),ncol=ncol(d)) #d is a binary matrix with 0.5's rounded down to 0's
    e_exc=(d)%*%t(!c)
    e_inc=(d)%*%t(d) #Need to use the version where 0.5 is cooerced to 0
    min_exc=matrix(mapply("min",e_exc,t(e_exc)),ncol=ncol(e_exc))
    out=matrix(mapply("min",min_exc,e_inc),ncol=ncol(e_exc))
    mean(out)
  }
  
  #Comparing the score to random shuffles at each locus
  dat=disagreement_score(filtered_muts$Genotype_shared_bin)
  nrand=100
  rand=sapply(1:nrand,function(i) {print(i);gt_shuffled<-apply(filtered_muts$Genotype_shared_bin,1,FUN = base::sample);disagreement_score(gt_shuffled)})
  
  sample_results=list(data=dat,rand=rand)
  
  all_results[[i]]<-sample_results
}

save(all_results,file="disagreement_scores")