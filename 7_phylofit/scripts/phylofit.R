#Script to run phylofit
#Run interactively on farm

#Author: Nick Williams (adapted by Emily Mitchell)

ID = "KX003"
Iteration = "KX003_5"
Age=81

setwd(paste0("/lustre/scratch117/casm/team154/em16/",ID,"/",Iteration,"/phylofit"))
#copy ultrametric tree and expanded clades file to this directory
#scp -r /Users/em16/Documents/PhD/Sequencing_results/DNA_seq/KX003/KX003_5/trees/vaf_cut/driver_info/KX003_5_expanded_clades_cut.csv em16@farm5-login:/lustre/scratch117/casm/team154/em16/KX003/KX003_5/phylofit
#scp -r /Users/em16/Documents/PhD/Sequencing_results/DNA_seq/KX003/KX003_5/trees/vaf_cut/tree_files/KX003_5_ultra.tree em16@farm5-login:/lustre/scratch117/casm/team154/em16/KX003/KX003_5/phylofit

library("rsimpop")
#library("kableExtra")
source("/lustre/scratch117/casm/team154/em16/R_scripts/functions/phylofit.R")

##Read in tree
tree=read.tree(paste0("/lustre/scratch117/casm/team154/em16/",ID,"/",Iteration,"/phylofit/",Iteration,"_ultra.tree"))
##Convert to years pcy (should I add 9months?)

tree$edge.length=tree$edge.length*Age
##Read in nodes..
dat=read.table(paste0("/lustre/scratch117/casm/team154/em16/",ID,"/",Iteration,"/phylofit/",Iteration,"_expanded_clades_cut.csv"),head=T,sep=",")
##Filter out clades with less than 5 samples.. Revisit this
dat=dat %>% filter(number_samples>3)
print(dat)

ud=dat %>% filter(grepl("^Clade",variant_ID))
##Super basic tree annotation: Color tree nodes..

pdf("drivers_tree.pdf", width=12,height=4)
tree=plot_tree(tree,b_do_not_plot = TRUE)
# blue for founder branch of established driver, red for unknown driver and black for other branches.
tree$coords$color=ifelse(tree$edge[,2] %in% dat$node,ifelse(tree$edge[,2] %in% ud$node,"blue","red"),"lightgrey")
tree=plot_tree(tree,cex.label =0,left.margin.prop = 0.1)
numcolonies=length(tree$tip.label)
# Configure tree annotation info
control=dat
control$cv=ifelse(grepl("^Clade",dat$variant_ID),"Expanded clade","Known Driver")
control$col=ifelse(grepl("^Clade",dat$variant_ID),"blue","red")
control$pch=15
control$label=control$variant_ID
control = list(vars=control,b.add.marker=TRUE,b.add.label=TRUE,cex.label=0.8)
for(node in dat$node){
  rsimpop:::add_specified_labels(list(details=data.frame(node=1)),
                                 tree = tree,
                                 node=node,
                                 control=control)
}
leg=with(unique(control$vars[,c("cv","col","pch")]),legend("topleft",legend=cv,col=col,pch=pch,pt.bg=col)$rect)
dev.off()

#Fit tree with and without final cell fraction

pdf("selection_coeff_plots.pdf", width=12, height=12)
par(layout = matrix(1:9, nrow=3, byrow=T))
out=lapply(1:dim(dat)[1],function(i) {
  #plot_tree(t1$ultratree,cex.label = 0,mar=par("mar"))
  nh=nodeHeights(tree)
  ht=nh[which(tree$edge[,2]==dat$node[i]),2]
  t1=drop.tip(tree,setdiff(tree$tip.label,get_samples_in_clade(dat$node[i],tree)))
  t1=add_outgroup(t1,num.shared.var = ht)
  t1$coords=NULL
  plot_tree(t1,cex.label = 0,mar=par("mar"))
  ## Fit without taking into account final aberrant cell fraction
  res1=fit_clade(tree,dat$node[i],nmutcolony = -1,nwtcolony =-1,maxt=1.5*Age)
  plot_res(res = res1,Age,b_add_extra = TRUE);title(sprintf('%s:No ACF Target',dat$variant_ID[i]))
  ## Fit taking into account final aberrant cell fraction
  res2=fit_clade(tree,dat$node[i],nmutcolony = dat$number_samples[i],nwtcolony = numcolonies-dat$number_samples[1],maxt=1.5*Age)
  plot_res(res = res2,Age,b_add_extra = TRUE);title(sprintf('%s:With ACF Target',dat$variant_ID[i]))
  list(res1=res1,res2=res2)
})
dev.off()



squan1=exp(t(sapply(out,function(x) quantile(x$res1$posterior$S,c(0.025,0.5,0.975)))))-1
squan2=exp(t(sapply(out,function(x) quantile(x$res2$posterior$S,c(0.025,0.5,0.975)))))-1
colnames(squan1)=c("S_lb","S","S_ub")
colnames(squan2)=c("S_lb","S","S_ub")
dat1=cbind(dat,squan1[,c("S","S_lb","S_ub")])
dat2=cbind(dat,squan2[,c("S","S_lb","S_ub")])

#Table of results with no acf target
write.table(dat1,file="res_no_acf_target.txt",sep="\t",col.names = TRUE,row.names = FALSE,quote=FALSE)
#kable(dat1, "html") %>% kable_styling("striped")

#Table of results with acf target
write.table(dat2,file="res_with_acf_target.txt",sep="\t",col.names = TRUE,row.names = FALSE,quote=FALSE)
#kable(dat2, "html") %>% kable_styling("striped")

#Copy to local driver
#scp -r em16@farm5-login:/lustre/scratch117/casm/team154/em16/KX003/KX003_5/phylofit/* /Users/em16/Documents/PhD/Sequencing_results/DNA_seq/KX003/KX003_5/phylofit/clade_size5
