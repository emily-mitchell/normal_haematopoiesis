#Author: Mike Spencer Chapman

library(stringr)
library(ape)
library(seqinr)
library(ggtree)
library(tidyr)
library(dplyr)
library(ggplot2)
library(gridExtra)
library(plotrix)
library(phangorn)
library(RColorBrewer)
library(Quartet)

#Define function required later in script for comparing phylogenies
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

get_RF_dist=function(tree1,tree2){
  comp<-comparePhylo(tree1,tree2)
  RF_dist=1-(length(comp$NODES[,1])/length(unique(tree$edge[,1])))
  return(RF_dist)
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

my_theme=theme_classic(base_family="Arial")+theme(text=element_text(size=7,family="Arial"),
                                                              axis.text=element_text(size=7,family="Arial"),
                                                              strip.text = element_text(size=7,family="Arial"),
                                                              legend.text = element_text(size=7,family="Arial"),
                                                              axis.title=element_text(size=7,family="Arial"),
                                                  axis.title.y=element_text(size=7,family="Arial"),
                                                  axis.title.x=element_text(size=7,family="Arial"))

my_working_directory="~/R_work/Phylogeny_of_foetal_haematopoiesis/"
treemut_dir="~/R_work/treemut" #Path for the cloned "treemut" directory
setwd(my_working_directory)

#Define the file paths for the data files
tree_file_path="Data/8pcw/Tree_8pcw.tree"
file_annot="Data/8pcw/Filtered_mut_set_annotated_8pcw"
bootstrapped_tree_stats_file_path="Data/8pcw/Bootstrapped_tree_stats"

R_function_files = list.files("R_functions",pattern=".R",full.names=TRUE)
sapply(R_function_files,source)
setwd(treemut_dir); source("treemut.R"); setwd(my_working_directory)

#Read count bootstraps
bootstrapped_trees_dir="~/Mounts/Lustre/fetal_HSC/fetal_8wks/tree_bootstraps/final_trees/"
all_bootstrapped_trees = "~/Mounts/Lustre/fetal_HSC/fetal_8wks/tree_bootstraps/all_trees"
bootstrap_tree_files=list.files(path=bootstrapped_trees_dir,full.names = T,pattern=".tree")
all_boots=read.tree(all_bootstrapped_trees)
all_boots<-di2multi(all_boots)

#Mutation bootstraps
bootstrapped_mutations_trees_dir="~/Mounts/Lustre/fetal_HSC/fetal_8wks/tree_bootstraps_mutations/final_trees/"
all_bootstrapped_mutations_trees = "~/Mounts/Lustre/fetal_HSC/fetal_8wks/tree_bootstraps_mutations/all_trees_8pcw_mutation_bootstraps"
bootstrap_mutations_tree_files=list.files(path=bootstrapped_mutations_trees_dir,full.names = T,pattern=".tree")
all_boots_mutations=read.tree(all_bootstrapped_mutations_trees)
all_boots_mutations<-di2multi(all_boots_mutations)

#Load up the main tree
tree=read.tree(tree_file_path)
ROOT=1+length(tree$tip.label)
tree_node_numbers=unique(tree$edge[,1]) #store list of the tree node numbers (not including the tips)

#Do comparisons for both sets of bootstraps - Quartet status
boots_comparisons_quartets=Quartet::QuartetStatus(all_boots,cf=tree)
boots_mutations_comparisons_quartets=Quartet::QuartetStatus(all_boots_mutations,cf=tree)

Quartet_divergence_plot<-cbind(SimilarityMetrics(boots_comparisons_quartets),type="8pcw: Sequencing read bootstraps")%>%
  rbind(cbind(SimilarityMetrics(boots_mutations_comparisons_quartets),type="8pcw: Mutation bootstraps"))%>%
  mutate(type=factor(type,levels=c("8pcw: Sequencing read bootstraps","8pcw: Mutation bootstraps")))%>%
  dplyr::select(QuartetDivergence,type)%>%
  gather(key="metric",value="value",-type)%>%
  ggplot(aes(x=value,fill=type))+
  geom_density(size=0.2)+
  facet_grid(rows=vars(type),scales="free_y")+
  #scale_x_continuous(limits=c(0.88,1))+
  my_theme + theme(legend.position = "none") +
  labs(x="Similarity score",y="Density")

boots_comparisons_splits=Quartet::SplitStatus(all_boots,cf=tree)
boots_mutations_comparisons_splits=Quartet::SplitStatus(all_boots_mutations,cf=tree)
RF_dist_plot<-data.frame(RF_dist=RobinsonFoulds(boots_mutations_comparisons_splits) / boots_mutations_comparisons_splits[, 'N'],type="8pcw: Mutation bootstraps")%>%
  rbind(data.frame(RF_dist=RobinsonFoulds(boots_comparisons_splits) / boots_comparisons_splits[, 'N'],type="read_count_bootstraps"))%>%
  mutate(type=factor(type,levels=c("8pcw: Sequencing read bootstraps","8pcw: Mutation bootstraps")))%>%
  ggplot(aes(x=RF_dist,fill=type))+
  geom_density(size=0.2)+
  facet_grid(rows=vars(type))+
  scale_x_continuous(limits=c(0,0.2))+
  my_theme + theme(legend.position = "none") +
  labs(x="Robinson-Foulds distance",y="Density")

RF_sim_plot<-cbind(RF_sim=RobinsonFoulds(boots_mutations_comparisons_splits,similarity = T)/boots_mutations_comparisons_splits[, 'N'],SimilarityMetrics(boots_mutations_comparisons_splits),type="8pcw: Mutation bootstraps")%>%
  rbind(cbind(RF_sim=RobinsonFoulds(boots_comparisons_splits,similarity = T)/boots_comparisons_splits[, 'N'],SimilarityMetrics(boots_comparisons_splits),type="8pcw: Sequencing read bootstraps"))%>%
  mutate(type=factor(type,levels=c("8pcw: Sequencing read bootstraps","8pcw: Mutation bootstraps")))%>%
  dplyr::select(RF_sim,type)%>%
  gather(key="metric",value="value",-type)%>%
  ggplot(aes(x=value,fill=type))+
  geom_density(size=0.2)+
  facet_grid(rows=vars(type),cols = vars(metric))+
  scale_x_continuous(limits=c(0.8,1))+
  my_theme + theme(legend.position = "none") +
  labs(x="Robinson-Foulds distance",y="Density")

QD_df<-cbind(SimilarityMetrics(boots_comparisons_quartets),type="8pcw: Sequencing read bootstraps")%>%
  rbind(cbind(SimilarityMetrics(boots_mutations_comparisons_quartets),type="8pcw: Mutation bootstraps"))%>%
  mutate(type=factor(type,levels=c("8pcw: Sequencing read bootstraps","8pcw: Mutation bootstraps")))%>%
  dplyr::select(QuartetDivergence,type)%>%
  gather(key="metric",value="value",-type)

RF_df<-cbind(RF_sim=RobinsonFoulds(boots_mutations_comparisons_splits,similarity = T)/boots_mutations_comparisons_splits[, 'N'],SimilarityMetrics(boots_mutations_comparisons_splits),type="8pcw: Mutation bootstraps")%>%
  rbind(cbind(RF_sim=RobinsonFoulds(boots_comparisons_splits,similarity = T)/boots_comparisons_splits[, 'N'],SimilarityMetrics(boots_comparisons_splits),type="8pcw: Sequencing read bootstraps"))%>%
  mutate(type=factor(type,levels=c("8pcw: Sequencing read bootstraps","8pcw: Mutation bootstraps")))%>%
  dplyr::select(RF_sim,type)%>%
  gather(key="metric",value="value",-type)

QD_and_RF_plot<-rbind(QD_df,RF_df)%>%
  filter(type=="8pcw: Sequencing read bootstraps")%>%
  ggplot(aes(x=value))+
  geom_density(size=0.2,fill=c("#E30613","#1D71B8")[1])+
  facet_grid(row = vars(metric),scales="free")+
  scale_x_continuous(limits=c(0.85,1))+
  #scale_y_continuous(limits=c(0,100))+
  my_theme + theme(legend.position = "none") +
  labs(x="Bootstrap tree similarity",y="Density")

comb_plots=arrangeGrob(Quartet_divergence_plot,RF_dist_plot,ncol=2,widths = c(2,1))
plot(comb_plots)
ggsave(comb_plots,filename = "~/Documents/Foetal_paper_revisions/Bootstrap_score_comparisons_8pcw.pdf",device=cairo_pdf,width = 7,height=3)
ggsave(QD_and_RF_plot,filename = "~/Documents/Foetal_paper_revisions/Bootstrap_readcount_metric_comparisons_8pcw.pdf",device=cairo_pdf,units="cm",width = 5,height=6)

#Create list of comparison stats to compare each of the bootstrap trees
boot_stats_list=list(nmuts=numeric(),
                     mean_mut_burden=numeric(),
                     nnodes_boot=numeric(),
                     shared_clades=list(),
                     new_clades=list())

#Fill up the stats list from the data
for(i in 1:length(all_boots)) {
  print(i)
  tree_boot=all_boots[[i]]
  comp=compare_nodes(tree,tree_boot) #Run the tree comparison function (hacked from the ape "comparePhylo" function)
  
  #Extract parameters of interest
  tree_boot_node_numbers=unique(tree_boot$edge[,1])
  tree_boot_lost_clades=tree_boot_node_numbers[which(!tree_boot_node_numbers%in%comp$tree_boot)]
  tree_boot_new_clades=lapply(tree_boot_lost_clades,function(node) {getTips(tree_boot,node)})
  
  #Store the output
  boot_stats_list$nmuts[i]<-sum(tree_boot$edge.length)
  boot_stats_list$mean_mut_burden[i]<-mean(get_mut_burden(tree_boot))
  boot_stats_list$nnodes_boot[i]<-length(tree_boot_node_numbers)
  boot_stats_list$shared_clades[[i]]<-tree_shared_clades
  boot_stats_list$tree_boot_new_clades[[i]]<-tree_boot_new_clades
}

#Create list of comparison stats to compare each of the bootstrap trees
boot_mutations_stats_list=list(nmuts=numeric(),
                     mean_mut_burden=numeric(),
                     nnodes_boot=numeric(),
                     shared_clades=list(),
                     new_clades=list())

#Fill up the stats list from the data
for(i in 1:length(all_boots_mutations)) {
  print(i)
  tree_boot=all_boots_mutations[[i]] #Load up the bootstrap tree; make multi-furcating
  comp=compare_nodes(tree,tree_boot) #Run the tree comparison function from ape "comparePhylo"
  
  #Extract parameters of interest
  tree_boot_node_numbers=unique(tree_boot$edge[,1])
  tree_boot_lost_clades=tree_boot_node_numbers[which(!tree_boot_node_numbers%in%comp$tree_boot)]
  tree_boot_new_clades=lapply(tree_boot_lost_clades,function(node) {getTips(tree_boot,node)})
  
  #Store the output
  boot_mutations_stats_list$nmuts[i]<-sum(tree_boot$edge.length)
  boot_mutations_stats_list$mean_mut_burden[i]<-mean(get_mut_burden(tree_boot))
  boot_mutations_stats_list$nnodes_boot[i]<-length(tree_boot_node_numbers)
  boot_mutations_stats_list$shared_clades[[i]]<-comp$tree
  boot_mutations_stats_list$tree_boot_new_clades[[i]]<-tree_boot_new_clades
}

save(boot_stats_list,boot_mutations_stats_list,file = bootstrapped_tree_stats_file_path)

load(bootstrapped_tree_stats_file_path)

#Save a sample of 10 of the bootstrap trees
pdf("~/Documents/Foetal_paper_revisions/Bootstrap_tree_comparisons_8pcw.pdf",width=16,height=5)
for(i in 1:10) {
  print(i)
  tree_boot=di2multi(read.tree(bootstrap_tree_files[i]))
  comparePhylo_and_plot(tree1=tree,tree2 = tree_boot,names=c("Original Tree (8pcw)","Bootstrap tree"))
}
dev.off()

retained_node_read_bootstraps=sapply(tree_node_numbers,function(node) {
  sum(unlist(lapply(boot_stats_list$shared_clades,function(x) node%in%as.numeric(x))))/length(boot_stats_list$nmuts)
})

retained_node_mutation_bootstraps=sapply(tree_node_numbers,function(node) {
  sum(unlist(lapply(boot_mutations_stats_list$shared_clades,function(x) node%in%as.numeric(x))))/length(boot_mutations_stats_list$nmuts)
})

RF_dist_read_bootstraps=sapply(boot_stats_list$shared_clades,function(x) 1-(length(x)/length(tree_node_numbers)))
RF_dist_mutation_bootstraps=sapply(boot_mutations_stats_list$shared_clades,function(x) 1-(length(x)/length(tree_node_numbers)))

p1<-data.frame(RF_dist=RF_dist_read_bootstraps,type="8pcw: Sequencing read bootstraps")%>%
  rbind(data.frame(RF_dist=RF_dist_mutation_bootstraps,type="8pcw: Mutation bootstraps"))%>%
  mutate(type=factor(type,levels=c("8pcw: Sequencing read bootstraps","8pcw: Mutation bootstraps")))%>%
  ggplot(aes(x=RF_dist,fill=type))+
  geom_density()+
  scale_x_continuous(limits=c(0,0.25))+
  facet_grid(rows=vars(type))+
  my_theme+theme(legend.position = "none")+
  labs(#title="Robinson-Foulds distance of bootstrapped trees",
    x="Robinson-Foulds distance",
    y="Density")

p2<-data.frame(number_of_nodes=boot_stats_list$nnodes_boot,type="8pcw: Sequencing read bootstraps")%>%
  rbind(data.frame(number_of_nodes=boot_mutations_stats_list$nnodes_boot,type="8pcw: Mutation bootstraps"))%>%
  mutate(type=factor(type,levels=c("8pcw: Sequencing read bootstraps","8pcw: Mutation bootstraps")))%>%
  ggplot(aes(x=number_of_nodes,fill=type))+
  geom_density()+
  geom_vline(xintercept=length(unique(tree$edge[,1])),col="red")+
  facet_grid(rows=vars(type))+
  my_theme+theme(legend.position = "none")+
  labs(#title = "Number of nodes in bootstrapped trees",
    x="Number of nodes",
    y="Density")

p3<-data.frame(number_of_mutations=boot_stats_list$nmuts)%>%
  ggplot(aes(x=number_of_mutations))+
  geom_density(fill="#9ECAE1")+
  geom_vline(xintercept=sum(tree$edge.length),col="red")+
  scale_x_continuous(limits=c(3600,4200))+
  my_theme+theme(legend.position = "none")+
  labs(#title = "Number of mutations in bootstrapped trees",
    x="Total mutations called",
    y="Density")

p4<-data.frame(mean_mutation_burden=boot_stats_list$mean_mut_burden)%>%
  ggplot(aes(x=mean_mutation_burden))+
  geom_density(fill="#9ECAE1")+
  geom_vline(xintercept=mean(get_mut_burden(tree)),col="red")+
  scale_x_continuous(breaks=seq(26,30,1),limits=c(25,30))+
  my_theme+
  labs(#title = "Mean mutation burden of bootstrapped trees",
    x="Mean mutation burden",
    y="Density")

comb_plot1=arrangeGrob(p1,p2,p3,p4,ncol=4)
plot(comb_plot1)
ggsave(comb_plot1,filename = "~/Documents/Foetal_paper_revisions/bootstrap_metrics_8pcw.pdf",device=cairo_pdf,width=7,height=2,units="in")


#Look specifically at the nodes at the top of the tree to test their robustness
first_three_division_nodes=c(279,280,314,315,317,318,380,402,403,447)
retained_node_read_bootstraps[tree_node_numbers%in%first_three_division_nodes]
retained_node_mutation_bootstraps[tree_node_numbers%in%first_three_division_nodes]

#KEY STATS FOR THE EARLY CLADES
#All but one of the clades in present in >95% of the bootstrapped trees
#The final clade, present in 78%, is that defined by one indel, which was not retained by filtering in >20% of bootstraps

data.frame(node=tree_node_numbers,retained=retained_node_read_bootstraps)%>%
  mutate(not_retained=1-retained)%>%
  mutate(node=reorder(node,not_retained))%>%
  mutate(early_node=ifelse(node%in%first_three_division_nodes,"yes","no"))%>%
  gather(key="Clade_status",value="proportion",-node,-early_node)%>%
  ggplot(aes(x=node,y=proportion,fill=Clade_status,col=factor(early_node),size=factor(early_node))) +
  geom_bar(stat="identity",width=0.85) +
  scale_fill_manual(values = c("#9ECAE1","#3182BD"))+
  scale_color_manual(values=c("black","red"),guide="none",)+
  scale_size_discrete(range=c(0.1,0.4),guide="none")+
  scale_y_continuous(breaks=seq(0,1,0.1))+
  my_theme+
  theme(axis.text.x = element_blank(),axis.ticks.x = element_blank())+
  labs(x="Original tree clades",
       y="Proportion of bootstraps where clade is retained",
       title = "Original tree clades retained in bootstrapped trees")

#KEY STATS
#40% of clades retained in 100% of bootstraps
#73% of clades retained in at least 95% of bootstraps
#83% of clades retained in at least 80% of bootstraps
#97.5% of clades retained in at least 50% of bootstraps

#These are the mpboot bootstrap support numbers
mpboot_df=data.frame(node=(length(tree$tip.label)+1):max(tree$edge),type="8pcw: MPBoot bootstrap approximation",retained=as.numeric(tree$node.label)/100)
boots_mutations_df=data.frame(node=tree_node_numbers,type="8pcw: Mutation bootstraps",retained=retained_node_mutation_bootstraps)

#Compare the results of read count bootstraps vs mutation bootstraps


p5.1<-data.frame(node=tree_node_numbers,type="8pcw: Sequencing read bootstraps",retained=retained_node_read_bootstraps)%>%
  mutate(node=reorder(node,-retained))%>%
  rbind(boots_mutations_df)%>%
  rbind(mpboot_df)%>%
  filter(node!=ROOT)%>%
  mutate(not_retained=1-retained)%>%
  mutate(early_node=ifelse(node%in%first_three_division_nodes,"yes","no"))%>%
  gather(key="Clade_status",value="proportion",-node,-early_node,-type)%>%
  mutate(type=factor(type,levels = c("8pcw: Sequencing read bootstraps","8pcw: Mutation bootstraps","8pcw: MPBoot bootstrap approximation")))%>%
  ggplot(aes(x=node,y=proportion,fill=Clade_status,col=factor(early_node),size=factor(early_node))) +
  geom_bar(stat="identity",width=0.85) +
  scale_fill_manual(values = c("#9ECAE1","#3182BD"))+
  scale_color_manual(values=c("black","red"),guide="none",)+
  scale_size_discrete(range=c(0.08,0.15),guide="none")+
  scale_y_continuous(breaks=seq(0,1,0.2))+
  facet_wrap(~type,nrow = 3)+
  my_theme+
  theme(axis.text.x = element_blank(),axis.ticks.x = element_blank(),legend.position = "none")+
  labs(x="Original tree clades",
       y="Proportion of bootstraps where clade is retained")

p5.2<-data.frame(node=tree_node_numbers,type="8pcw: Sequencing read bootstraps",retained=retained_node_read_bootstraps)%>%
  mutate(node=reorder(node,-retained))%>%
  rbind(boots_mutations_df)%>%
  filter(node!=ROOT)%>%
  mutate(not_retained=1-retained)%>%
  mutate(early_node=ifelse(node%in%first_three_division_nodes,"yes","no"))%>%
  gather(key="Clade_status",value="proportion",-node,-early_node,-type)%>%
  mutate(type=factor(type,levels = c("8pcw: Sequencing read bootstraps","8pcw: Mutation bootstraps")))%>%
  ggplot(aes(x=node,y=proportion,fill=Clade_status,col=factor(early_node),size=factor(early_node))) +
  geom_bar(stat="identity",width=0.85) +
  scale_fill_manual(values = c("#9ECAE1","#3182BD"))+
  scale_color_manual(values=c("black","red"),guide="none",)+
  scale_size_discrete(range=c(0.1,0.3),guide="none")+
  scale_y_continuous(breaks=seq(0,1,0.2))+
  facet_wrap(~type,nrow = 3)+
  my_theme+
  theme(axis.text.x = element_blank(),axis.ticks.x = element_blank(),legend.position = "none")+
  labs(x="Original tree clades",
       y="Proportion of bootstraps where clade is retained",
       title = "Original tree clades retained in bootstrapped trees")

#Highlight the lower confidence nodes - low confidence by all bootstrap methods
low_confidence_nodes=cbind(data.frame(node=tree_node_numbers,
                 read_bootstraps=retained_node_read_bootstraps,
                 mutation_bootstraps=retained_node_mutation_bootstraps,
                 mpboot_bootstrap_approximation=as.numeric(tree$node.label)/100))%>%
  filter(read_bootstraps<0.9 &mutation_bootstraps<0.9 & mpboot_bootstrap_approximation<0.9)
cairo_pdf(filename = "~/Documents/Foetal_paper_revisions/Lower_confidence_split_8pcw.pdf",width = 15,height=8)
tree=plot_tree(tree,cex.label=0)
sapply(1:nrow(low_confidence_nodes), function(i) {
  node=low_confidence_nodes$node[i]
  confidence=mean(low_confidence_nodes$read_bootstraps[i],low_confidence_nodes$mutation_bootstraps[i],low_confidence_nodes$mpboot_bootstrap_approximation[i])
  info=get_edge_info(tree,details=NULL,node)
  col.func=colorRampPalette(colors=c("red","blue"))
  cols=col.func(100)
  draw.circle(x=info$x,y=info$yb,radius=1,col = cols[100*round(confidence,digits=2)],lwd=0.1)
})
dev.off()

p6<-data.frame(node=tree_node_numbers,retained=retained_node_read_bootstraps)%>%
  filter(node!=ROOT)%>%
  mutate(clade_size=sapply(node,function(node) length(getTips(tree,node))),
         nmuts=sapply(node,function(node) tree$edge.length[tree$edge[,2]==node]))%>%
  mutate(nchances=clade_size*nmuts)%>%
  ggplot(aes(x=nchances,y=retained))+
  geom_jitter(col="#3182BD",alpha=0.5,size=2)+
  scale_x_log10(breaks=c(2,4,6,8,16,32,64,128,254))+
  scale_y_continuous(breaks=seq(0,1,0.1))+
  my_theme+
  labs(x="Number of samples in clade x Number of mutations defining clade",
       y="Proportion of bootstraps where clade is retained",
       title = "Clade robustness in bootstrapping")

#Extract information on the new clades - find all unique new clades, and define their frequencies in the bootstraps
all_new_clades=unique(unlist(lapply(boot_stats_list$tree_boot_new_clades,function(list) {sapply(list,paste,collapse=",")})))
new_clades_by_bootstrap<-lapply(boot_stats_list$tree_boot_new_clades,function(list) {sapply(list,paste,collapse=",")})
new_clade_freqs=sapply(all_new_clades,function(clade) {sum(unlist(lapply(new_clades_by_bootstrap,function(x) clade%in%x)))})

#Plot the frequency of these new clades
p7<-data.frame(clade=1:length(new_clade_freqs),clade_freq=new_clade_freqs/length(boot_stats_list$nnodes_boot))%>%
  filter(clade_freq>0.01)%>% #Only include clades found in >1% of the bootstraps (there is a long tail of rare clades)
  mutate(clade=reorder(clade,-clade_freq))%>%
  ggplot(aes(x=clade,y=clade_freq))+
  geom_bar(stat="identity",col="black",size=0.1,fill="#9ECAE1")+
  scale_y_continuous(limits=c(0,1),breaks=seq(0,1,0.1))+
  my_theme+
  theme(axis.text.x = element_blank(),axis.ticks.x = element_blank())+
  labs(x="Novel clades",
       y="Proportion of bootstraps containing novel clade",
       title = "Novel clades found in bootstrapped trees")

comb_p6_p7=arrangeGrob(p7,p6,ncol=2)
comb_plot2<-arrangeGrob(p5.1,comb_p6_p7,nrow=2,heights = c(2,1))
plot(comb_plot2)
ggsave(comb_plot2,filename = "~/Documents/Foetal_paper_revisions/Bootstrap_trees_clade_robustness_8pcw.pdf",width = 10,height=10,units = "cm",device=cairo_pdf)
ggsave(p5.1,filename = "~/Documents/Foetal_paper_revisions/Bootstrap_3methods_trees_clade_robustness_8pcw.pdf",width = 9,height=8.5,units="cm",device=cairo_pdf)
