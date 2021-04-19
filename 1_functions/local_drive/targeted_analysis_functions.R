get_sample_mutations = function(sample,tree,details,vcf_file=FALSE) {
  terminal_node <- which(tree$tip.label == sample)
  all_nodes <- get_ancestral_nodes(terminal_node,tree$edge)
  mutations <- details[details$node %in% all_nodes,]
  if(vcf_file==TRUE) {
    vcf_file = create_vcf_files(mat=mutations)
    return(vcf_file)
  } else {
    return(mutations)
  }
}

validate_colony_muts=function(colony,tree,details,NV,NR,validatable_depth_cutoff=8,pval_cutoff=0.05) {
  tree$tip.label=gsub("_hum","",tree$tip.label)
  colony_muts = get_sample_mutations(colony,tree,details = details)$mut_ref
  n_colony = length(colony_muts)
  NR_colony=NR[colony_muts,colony]
  NR_colony[NR_colony==0] <- 1
  NV_colony=NV[colony_muts,colony]
  n_colony_validatable = sum(NR_colony >= 8)
  p_colony_validatable = n_colony_validatable/n_colony
  pvals = sapply(1:length(NV_colony), function(i) {binom.test(NV_colony[i],NR_colony[i],p=ifelse(grepl("X",colony_muts[i])|grepl("Y",colony_muts[i]),0.95,0.5),alternative="less")$p.value})
  pvals_v = pvals[NR_colony>=8]
  n_colony_validated = sum(pvals_v > 0.02)
  p_colony_validated = n_colony_validated/n_colony_validatable
  
  private_muts=details_targ_full[details$node==which(tree$tip.label==colony),"mut_ref"]
  n_private = length(private_muts)
  if(n_private>0){
    NR_private = NR[private_muts,colony]
    NR_private[NR_private==0] <-1
    NV_private = NV[private_muts,colony]
    n_private_validatable = sum(NR_private >= validatable_depth_cutoff)
    p_private_validatable = n_private_validatable/n_private
    pvals = sapply(1:length(NV_private), function(i) {binom.test(NV_private[i],NR_private[i],p=ifelse(grepl("X",private_muts[i])|grepl("Y",private_muts[i]),0.95,0.5),alternative="less")$p.value})
    pvals_v = pvals[NR_private>=validatable_depth_cutoff]
    n_private_validated = sum(pvals_v > pval_cutoff)
    p_private_validated = n_private_validated/n_private_validatable
    n_likely_subclonal=sum(pvals<pval_cutoff & NR_private>=validatable_depth_cutoff & NV_private>=3) #Good number of reads to support variant being real, but not consistent with true somatic probability distribution
  } else{
    n_private_validatable=n_private_validated=n_likely_subclonal=0
    p_private_validatable=p_private_validated=NA
  }
  colony_out=data.frame(sample=colony,
                        mean_cov=mean(NR[,colony]),
                        on_existing_phylogeny=colony%in%tree$tip.label,
                        n_colony=n_colony,
                        n_colony_validatable=n_colony_validatable,
                        p_colony_validatable=p_colony_validatable,
                        n_colony_validated=n_colony_validated,
                        p_colony_validated = p_colony_validated,
                        n_private = n_private,
                        n_private_validatable=n_private_validatable,
                        p_private_validatable=p_private_validatable,
                        n_private_validated=n_private_validated,
                        p_private_validated=p_private_validated,
                        n_likely_subclonal=n_likely_subclonal)
  return(colony_out)
}

validate_mutation = function(mutation,tree,details,NV,NR,pval_cutoff=0.05,vaf_cutoff=0.3,depth_cutoff=8,counts_only=FALSE) {
  require(ape)
  if(mutation %in% rownames(NV)){
    mut_node=details$node[details$mut_ref==mutation]
    if(mut_node>length(tree$tip.label)){
      samples=gsub(pattern="_hum",replacement="",extract.clade(tree,mut_node)$tip.label)
      if(any(samples %in% colnames(NV))) {
        NV_mut=sum(NV[mutation,samples[samples %in% colnames(NV)]])
        NR_mut=sum(NR[mutation,samples[samples %in% colnames(NV)]])
        if(counts_only) {stop(return(c(NV_mut,NR_mut)))}
      }else{
        validation_result <- "No targeted samples expected to have mutation"
        stop(return(validation_result))
      }
    }else{
      samples=gsub(pattern="_hum",replacement="",tree$tip.label[mut_node])
      if(any(samples %in% colnames(NV))){
        NV_mut=NV[mutation,samples]
        NR_mut=NR[mutation,samples]
        if(counts_only) {stop(return(c(NV_mut,NR_mut)))}
      }else{
        validation_result<-"No targeted samples expected to have mutation"
        stop(return(validation_result))
      }
    }
    if(NR_mut>=depth_cutoff) {
      pval=binom.test(NV_mut,NR_mut,p=ifelse(grepl("X",mutation)|grepl("Y",mutation),0.95,0.5),alternative="less")$p.value
      vaf=sum(NV_mut)/sum(NR_mut) #Use alternative vaf cutoff for shared mutations.  Shared mutations have higher depth, and even a small bias in capture/ sequencing rates can result in failing the binomial test.
      if(pval>pval_cutoff|(vaf>=vaf_cutoff&length(samples>1)&sum(NV_mut)>=10)){validation_result <- "PASS"} else {validation_result <- "FAIL"}
    }else{
      validation_result <- "Inadequate depth"
    }
  }else{
    validation_result<-"Not in bait set"
  }
  return(validation_result)
}

view_validation_vaf_plot=function(sample, tree, details,NV,NR,depth_cutoff=0,include_muts="all") {
  tree$tip.label=gsub("_hum","",tree$tip.label)
  if(include_muts=="all") {
    mutations=get_sample_mutations(sample,tree,details)$mut_ref
  }else if(include_muts=="private"){
    mutations=details[details$node==which(tree$tip.label==sample),"mut_ref"]
  }
  mutations<-mutations[mutations%in%rownames(NV)] #Include only mutations in the bait set
  mutations<-mutations[NR[mutations,sample]>=depth_cutoff]
  if(length(mutations)<2){stop(return("<2 mutations meeting criteria - indequate to create plot\n"))}
  n_fail=sum(details$validation_results[details$mut_ref %in% mutations] == "FAIL")
  n_pass=sum(details$validation_results[details$mut_ref %in% mutations] == "PASS")
  sample_mean_depth = mean(NR[mutations,sample])
  dens <- density((NV/NR)[mutations,sample])
  plot(dens, main = sample,xlim=c(-0.05,1.05))
  abline(v = dens$x[which.max(dens$y)])
  text(0.7,
       max(dens$y) - 0.2,
       paste("Peak VAF dens=",
             round(dens$x[which.max(dens$y)], digits = 2),
             "\nDepth cutoff for validation=",
             depth_cutoff,
             "\nMean coverage of included mutations=",
             round(sample_mean_depth,digits =2),
             "\nMutations included=",
             length(mutations),
             "\nFail variants:",
             n_fail,
             "\nPass variants:",
             n_pass), col = "red", cex = 0.7)
}

#These are functions required to calculate the true daughter nodes from a polytomous tree that is structured as a dichotomous tree
get_level_daughters = function(node,tree){
  ancestral_node_height=nodeHeights(tree)[tree$edge[,2]==node,2]
  all_daughters=get_all_node_children(node,tree)
  all_daughters_ordered=tree$edge[tree$edge[,2]%in%all_daughters,2]
  all_daughters_heights=nodeHeights(tree)[tree$edge[,2]%in%all_daughters,2]
  return(all_daughters_ordered[all_daughters_heights==ancestral_node_height])
}
get_direct_daughters=function(ancestral_node,tree){
  level_ancestors=c(ancestral_node,get_level_daughters(ancestral_node,tree_targ))
  all_daughters=tree$edge[tree$edge[,1]%in%level_ancestors,2]
  true_daughters=all_daughters[!all_daughters%in%level_ancestors]
  return(true_daughters)
}

#Print the "lineage loss stats" (how much the daughter node cell fractions add up to the parent node cell fraction)
print_lineage_loss_stats = function(ancestral_node,sample,tree,details,matrices=list()) {
  ancestral_node_frac = get_node_cell_frac(ancestral_node,sample,tree=tree,details=details,matrices=matrices)
  daughters=get_direct_daughters(ancestral_node,tree=tree_targ)
  daughter_node_fracs = rep(0,length(daughters))
  for(i in 1:length(daughters)) {
    cell_frac=get_node_cell_frac(daughters[i],sample=sample,details = details,tree=tree,matrices = matrices)
    daughter_node_fracs[i] <- cell_frac
  }
  names(daughter_node_fracs) <- daughters
  print("Parent node cell fraction:")
  print(round(ancestral_node_frac,digits=3))
  print("Individual daughter node cell fractions:")
  print(round(daughter_node_fracs,digits=3))
  total_daughter_frac=sum(daughter_node_fracs)
  print("Summed daughter node cell fraction is:")
  print(round(total_daughter_frac,digits=3))
  print("Daughters/ Parents ratio is:")
  print(round(as.numeric(total_daughter_frac/ancestral_node_frac),digits=3))
}

#Function to calculate the cell fraction in a bulk sample for the mutations leading to a specific node
get_node_cell_frac = function(node,sample,tree,details,matrices){
  node_muts=details$mut_ref[details$node==node]
  if(length(node_muts)==0){
    return(NA)
  } else {
    auto_muts=node_muts[!grepl("X",node_muts)&!grepl("Y",node_muts)]
    xy_muts=node_muts[grepl("X",node_muts)|grepl("Y",node_muts)]
    auto_cell_NV = sum(matrices$mtr[auto_muts,sample])
    auto_cell_NR = sum(matrices$dep[auto_muts,sample])
    xy_cell_NV = sum(matrices$mtr[xy_muts,sample])
    xy_cell_NR = sum(matrices$dep[xy_muts,sample])
    if(auto_cell_NR+xy_cell_NR == 0) {
      return(NA)
    } else {
      cell_frac_CI = binom.test(auto_cell_NV+xy_cell_NV,auto_cell_NR+2*xy_cell_NR)$conf.int[1:2] * 2
      overall_cell_frac=sum(c(auto_cell_NV,xy_cell_NV))/sum(c((auto_cell_NR/2),xy_cell_NR))
      attr(overall_cell_frac,"conf.int") <- cell_frac_CI
      return(overall_cell_frac)
    }
  }
}

#PLOT NODE LABELS FUNCTIONS (to add to Nick's functions)
#The "plot_node_cell_frac" function
plot_node_cell_frac = function(node,sample,tree,details,matrices,cex=0.6) {
  node_cell_frac=get_node_cell_frac(node,sample,tree,details,matrices)
  info=get_edge_info(tree,details,node)
  if(!is.na(node_cell_frac) & node_cell_frac>0.002) {
    text(info$x,info$yb,round(node_cell_frac,digits=3),cex = cex,col="black",font=2)
  }
}
#The "plot_node_number" function
plot_node_number = function(tree,details,matrices,node,cex=0.4) {
  info=get_edge_info(tree,details,node)
  text(info$x,info$yb,node,cex = cex,col="black",font=2)
}
#The add_annotation_targeted is for functions that require a SAMPLE name to list the bulk sample that is being overlayed on tree (and who's read counts are in the NV/ NR matrices)
add_annotation_targeted=function(sample,tree,details,matrices,annot_function,plot_sample_name=TRUE,...){
  N=dim(tree$edge)[1]
  lapply(1:N,function(i) annot_function(tree$edge[i,2],sample,tree,details,matrices,...))
  if(plot_sample_name) {text(20,2,sample)}
}

#Colour the mutations based on a categorial variable - in the first case a validation result
add_categorical_col=function(tree, ##<< enhanced phylo returned from plot_tree
                             details,##<< dataframe with summary details of mutations mapped to tree together with EDGE_IDX - index of mutation in edges matrix
                             matrices,##<< a list of matrices parallel to details with columns named by tip labels
                             node,
                             var_field,
                             annot=list(),
                             b.add.line=TRUE,
                             ...){
  ##Get all the detail about the edge coords + idx in detail
  info=get_edge_info(tree,details,node)
  if(length(info$idx.in.details) > 0) {
    bdat=details[[var_field]][info$idx]
    bdat_filt=bdat[bdat %in% names(annot)]
    if(is.null(bdat)){
      stop("Error in provided bfield (does it exist and is it numeric?)")
    }
    bdat_filt = sort(bdat_filt,decreasing=TRUE)
    if(length(bdat_filt)>0){
      edge_length=length(bdat_filt)
      ##Could add in a third category NA
      #missing=sum(is.na(bdat))
      if(b.add.line){
        y0_next = info$yt
        for(i in 1:edge_length) {
          arrows(y0=y0_next,y1=(y0_next - 1),x0=info$x,length = 0,col=annot[[bdat_filt[i]]],lend=1,...)
          y0_next = y0_next - 1
        }
      }
    }
  }
}

#Plot just the tip labels (or a point) for an individual sample
plot_sample_tip_label = function(tree,details,sample) {
  node=which(tree$tip.label==sample)
  info=get_edge_info(tree,details,node)
  text(sample,x=info$x,y=info$yb,cex = 0.7)
}

plot_sample_tip_point = function(tree,details,sample) {
  node=which(tree$tip.label==sample)
  info=get_edge_info(tree,details,node)
  points(x=info$x,y=info$yb,type="p",pch=21,bg="black",col="black")
}

plot_d_or_r_tip_point = function(sample,tree,details,donor_ID,recip_ID) {
  node=which(tree$tip.label==sample)
  info=get_edge_info(tree,details,node)
  tip_col=ifelse(grepl(donor_ID,sample),"dark green","red")
  points(x=info$x,y=info$yb,type="p",pch=20,bg=tip_col,col=tip_col)
}

#Nick's "get_ancestral_nodes" function. Needed for the validation functions above
get_ancestral_nodes= function(node,edge,exclude_root=TRUE){
  idx=which(edge[,2]==node)
  parents=node ##Include the node
  while(length(idx)>0){
    if(length(idx)>1){
      stop("multiple parents!")
    }
    parent=edge[idx,1]
    parents=c(parents,parent)
    #This finds the parent of the current parent - thus navigating up to the root.
    idx=which(edge[,2]==parent)
  }
  if(exclude_root){
    parents[-length(parents)] ##The last node is the root.
  }else{
    parents
  }
}

#Function to extract the node counts from mutations on a given branch, for a given sample
#Divides counts into auto & XY mutations, and returns as a df
get_node_read_counts=function(node,sample,tree,details,matrices,exclude_mut_indexes=NULL) {
  info=get_edge_info(tree.multi,details,node)
  
  #Exclude mutations from a branch if they match those given
  if(!is.null(exclude_mut_indexes)) {
    info$idx.in.details<-info$idx.in.details[!info$idx.in.details%in%exclude_mut_indexes]
  }
  
  #Separate out the indexes of the autosomal and XY chromosomes
  info$idx.in.details.AUTO=info$idx.in.details[!grepl("X",details$mut_ref[info$idx.in.details])&!grepl("Y",details$mut_ref[info$idx.in.details])]
  info$idx.in.details.XY=info$idx.in.details[grepl("X",details$mut_ref[info$idx.in.details])&!grepl("Y",details$mut_ref[info$idx.in.details])]
  
  #Store the read counts of each in the df
  NV_auto=sum(matrices$NV[info$idx.in.details.AUTO,sample])
  NR_auto=sum(matrices$NR[info$idx.in.details.AUTO,sample])
  NV_xy=sum(matrices$NV[info$idx.in.details.XY,sample])
  NR_xy=sum(matrices$NR[info$idx.in.details.XY,sample])
  
  #Return as a vector
  #return(c(NV_auto,NR_auto,NV_xy,NR_xy))
  
  #Or as a df
  df=data.frame(node=node,NV_auto=NV_auto,NR_auto=NR_auto,NV_xy=NV_xy,NR_xy=NR_xy)
  return(df)
}

#Function to bootstrap the counts (of format above), to return cell fraction bootstrapped estimates
bootstrap_counts = function(node_counts,boot_straps=1000) {
  NR_auto=node_counts$NR_auto;NV_auto=node_counts$NV_auto
  NR_xy=node_counts$NR_xy;NV_xy=node_counts$NV_xy
  if((NR_auto+NR_xy)==0) {
    cell_fracs<-rep(0,boot_straps)
  } else if(NR_xy==0){
    cell_fracs=(rbinom(boot_straps,NR_auto,prob = NV_auto/NR_auto))/(NR_auto/2)
  } else if(NR_auto==0) {
    cell_fracs=(rbinom(boot_straps,NR_xy,prob=NV_xy/NR_xy))/NR_xy
  } else {
    cell_fracs=(rbinom(boot_straps,NR_auto,prob = NV_auto/NR_auto) +
                  rbinom(boot_straps,NR_xy,prob=NV_xy/NR_xy))/sum(c((NR_auto/2),NR_xy))
  }
}

#Function to check the distribution for a given sample
check_branch_distribution=function(sample,node,tree,details,matrices,return_counts=FALSE) {
  #Get node info
  info=get_edge_info(tree,details,node)
  
  #Separate out the indexes of the autosomal and XY chromosomes
  info$idx.in.details.AUTO=info$idx.in.details[!grepl("X",details$mut_ref[info$idx.in.details])&!grepl("Y",details$mut_ref[info$idx.in.details])]
  info$idx.in.details.XY=info$idx.in.details[grepl("X",details$mut_ref[info$idx.in.details])&!grepl("Y",details$mut_ref[info$idx.in.details])]
  
  if(length(info$idx.in.details.AUTO)==0) {
    x_auto<-NULL
    n_auto<-NULL
  } else {
    x_auto=matrices$NV[info$idx.in.details.AUTO,sample]
    n_auto=matrices$NR[info$idx.in.details.AUTO,sample]
  }
  if(length(info$idx.in.details.XY)==0) {
    x_xy<-NULL
    n_xy<-NULL
  } else {
    x_xy=matrices$NV[info$idx.in.details.XY,sample]
    n_xy=matrices$NR[info$idx.in.details.XY,sample]
  }
  x=c(x_auto,x_xy)
  n=c(n_auto,2*n_xy)
  muts=c(details$mut_ref[info$idx.in.details.AUTO],details$mut_ref[info$idx.in.details.XY])
  idx=c(info$idx.in.details.AUTO,info$idx.in.details.XY)
  df=data.frame(mut_ref=muts,idx.in.details=idx,NV=x,NR=n)
  n[n==0]<-1
  if(length(x)>1) {
    pval=prop.test(x,n)$p.value
  } else{
    pval<-NA
  }
  names(pval)=sample
  if(return_counts){
    return(df)
  } else {
    return(pval)
  }
}

#Function to look through the muts in a branch and define counts that come from the earliest (or not earliest) binomial distribution
find_early_muts_from_branch=function(sample,node,tree,details,matrices,return_late_muts=FALSE,cols=c("red","orange","blue","purple")) {
  print(node)
  df_node=check_branch_distribution(sample,node,tree,details,matrices,return_counts = TRUE)
  
  #Perform binomial mixture model to work out which mutations are the early ones on each branch
  hist(df_node$NV/df_node$NR,col='gray',freq=F,xlab="VAF",main=paste("Node VAFs for node",node))
  lines(density(df_node$NV/df_node$NR),lwd=2,lty='dashed')
  
  if(nrow(df_node)==2) {
    max_p = which.max(df_node$NV/df_node$NR)
    early_idx=c(1:length(df_node$NV)%in%max_p)
  } else {
    res = binom_mix(df_node$NV,df_node$NR,nrange=2:min(3,length(df_node$NV)))
    for (i in 1:res$n){
      meancov = round(mean(df_node$NR))
      lines(x=(0:meancov)/meancov,
            y=meancov*res$prop[i]*dbinom(0:meancov,meancov,prob=res$p[i]),
            type="l",col=cols[i],lwd=2)
    }
    max_p=which.max(res$p)
    early_idx=res$Which_cluster==max_p
  }
  
  if(return_late_muts) {
    max_p_idx=df_node$idx.in.details[!early_idx]
  } else {
    max_p_idx=df_node$idx.in.details[early_idx]
  }
  return(max_p_idx)
}

node_lineage_loss = function(node,sample,tree,details,matrices,boot_straps,CI=0.95,display_vafs=FALSE,return_ancestral_cell_frac=FALSE) { #tree must be a multifurcating tree
  sample=gsub("_comb","",sample)
  ancestral_node_counts=get_node_read_counts(node,sample,tree,details,matrices)
  ancestral_cell_fracs = bootstrap_counts(ancestral_node_counts,boot_straps = boot_straps)
  if(display_vafs) {hist(ancestral_cell_fracs,main="Ancestral node cell fraction distribution");abline(v=median(ancestral_cell_fracs))}
  daughters=tree$edge[tree$edge[,1]==node,2]   #Find the daughter nodes - this is why must be multifurcating tree
  #Get counts from each individual node by apply lapply; then Reduce to a df
  daughter_node_counts=Reduce(rbind,lapply(daughters,get_node_read_counts,sample=sample,tree=tree,details=details,matrices=matrices))
  #Now bootstrap across the rows of the df to get df of bootstrapped cell fractions
  daughter_cell_fracs=Reduce(rbind,lapply(1:nrow(daughter_node_counts), function(i) {bootstrap_counts(daughter_node_counts[i,])}))
  
  if(display_vafs) {apply(daughter_cell_fracs,1,hist,main="Daughter cell fraction distributions")}
  
  #default value to output
  lineages_captured=colSums(daughter_cell_fracs)/ancestral_cell_fracs
  lineages_captured[lineages_captured>1]<-1
  if(display_vafs) {hist(lineages_captured,main = "Proportion of daughter lineages captured in phylogeny")}
  
  #Alternative stat for "cell fraction lost at node" - i.e. absolute values, not relative
  # cell_frac_lost=ancestral_cell_fracs-colSums(daughter_cell_fracs)
  # cell_frac_lost[cell_frac_lost<0]<-0
  # hist(cell_frac_lost)
  
  if(return_ancestral_cell_frac) {
    ancestral_cell_fracs[ancestral_cell_fracs>1] <- 1
    df=data.frame(median=median(ancestral_cell_fracs,na.rm=TRUE),lower_CI=quantile(ancestral_cell_fracs,(1-CI)/2,na.rm=TRUE),upper_CI=quantile(ancestral_cell_fracs,1-(1-CI)/2,na.rm=TRUE))
  } else {
    df=data.frame(median=median(lineages_captured,na.rm=TRUE),lower_CI=quantile(lineages_captured,(1-CI)/2,na.rm=TRUE),upper_CI=quantile(lineages_captured,1-(1-CI)/2,na.rm=TRUE))
  }
  return(df)
}

plotDonut=function(x,y,median=NA,radius,col,prop,llwd=0.5,border=NA,width=NA,plotPie=FALSE){
  lims=par("usr")
  as=dev.size()
  asr=as[1]/as[2]
  yscale=asr*(lims[4]-lims[3])/(lims[2]-lims[1])
  prop=prop/sum(prop)
  cutpoint=c(0,cumsum(prop)*2*pi)
  
  N=2*pi/0.05
  n=ceiling(N*diff(cutpoint)/(2*pi))
  d=diff(cutpoint)/n
  if(!any(is.na(cutpoint))&!any(is.nan(cutpoint))) {
    if(length(prop)>1){
      for(i in 2:length(cutpoint)){
        if(cutpoint[i-1]==cutpoint[i]) {next}
        polygon(x+c(radius*sin(seq(cutpoint[i-1],cutpoint[i],d[i-1])),0),
                y+yscale*c(radius*cos(seq(cutpoint[i-1],cutpoint[i],d[i-1])),0),
                border=border,col=col[i-1],lwd=llwd)
      }
    }else{
      i=2
      polygon(x+c(radius*sin(seq(cutpoint[i-1],cutpoint[i],d[i-1]))),
              y+yscale*c(radius*cos(seq(cutpoint[i-1],cutpoint[i],d[i-1]))),
              border=border,col=col[1],lwd=llwd)
    }
    yscale
    #Cut out the "donut" from the centre
    if(!plotPie) {
      polygon(x+c((radius/3)*cos(seq(from=0,to=2*pi,by=(2*pi/N)))),
              y+yscale*c((radius/3)*sin(seq(0,2*pi,by=(2*pi/N)))),
              border=border,col="white",lwd=llwd) 
    }
    #Now to draw the line for the median
    if(median>0 & !is.na(median)) {
      lines(x+c((0.66*radius)*sin(seq(from=0,to=2*median*pi,by=(2*pi/N)))),
            y+yscale*c((0.66*radius)*cos(seq(0,2*median*pi,by=(2*pi/N)))),
            col="black",lwd=3,lend=1)
    }
  }
}

###THE MAIN PLOTTING FUNCTION FOR THE TREES##
#Not entirely self-sufficient.  Needs various appropriately labelled objects in the environment

generate_targ_seq_plots=function(samples,
                                 tree,
                                 details_targ,
                                 matrices,
                                 post_prob_type=c("raw","clean"),
                                 info_type=c("post.prob","cell_frac","log_cell_frac"),
                                 prob_threshold_to_include=0.5, #Probability threshold from the post.prob matrix for plotting
                                 plot_cell_frac=TRUE,
                                 plot_donut=TRUE,
                                 donut_info="cell_frac", #other option is "lineages_lost"
                                 CI=0.8,  #Confidence intervals on the pie chart, default = 80% CI
                                 radius=3.5,  #Radius of the pie charts on the plot
                                 scale_muts_to_branch=FALSE) {
  require(plotrix)
  if(post_prob_type=="raw") {
    post.prob.mat <-post.prob[details_targ$mut_ref,]
  } else if(post_prob_type=="clean"){
    post.prob.mat <-clean.post.prob[details_targ$mut_ref,]
  }
  if(info_type=="post.prob") {
    post.prob.mat[post.prob.mat<0.05] <- 0
    details_targ_full=cbind(details_targ,post.prob.mat)
  } else {
    cell_frac=calculate_cell_frac(matrices$NV,matrices$NR)
    cell_frac_present<-cbind(cell_frac,cell_frac[,gsub("_comb","",colnames(post.prob.mat)[grep("_comb",colnames(post.prob.mat))])]) #To double up the "_comb" results, so that matched the post.prob.mat
    cell_frac_present=cell_frac_present[details_targ$mut_ref,]
    colnames(cell_frac_present)<- colnames(post.prob)
    cell_frac_present[post.prob.mat<prob_threshold_to_include]<-0
    if(info_type=="cell_frac") {
      details_targ_full=cbind(details_targ,cell_frac_present)
    } else if(info_type=="log_cell_frac") {
      log_cell_frac_present=cell_frac_present
      log_cell_frac_present[cell_frac_present != 0] <- log(log_cell_frac_present[cell_frac_present != 0])#change all the non 0 vafs to the log of their vaf
      log_cell_frac_present_scaled=log_cell_frac_present
      scale_range=c(0.01,1)
      log_cell_frac_present_scaled[cell_frac_present!=0] = plotrix::rescale(log_cell_frac_present[cell_frac_present!=0],newrange = scale_range) #scale these figures between 0 and 1
      details_targ_full=cbind(details_targ,log_cell_frac_present_scaled)
    }
  }
  lims=par("usr")
  sapply(samples, function(sample) {
    tree=plot_tree(tree, cex.label = 0,lwd=0.5,plot_axis = TRUE,default_edge_color="lightgrey")
    lims=par("usr")
    if(grepl("PD",sample)) {
      text(0,lims[4]-1,paste0(lcm_smry$Tissue[lcm_smry$Sample_ID==sample]," (",sample,"): Mean depth is ",round(mean(NR[,gsub("_comb","",sample)]),digits = 2)),cex=1,pos=4)
    } else {
      text(0,lims[4]-0.5,paste0(sample,": Mean depth is ",round(mean(NR[,gsub("_comb","",sample)]),digits = 2)),cex=1,pos=4)
    }
    if(grepl("PD",sample)) {
      text(0,lims[4]-3,paste0("Library concentration was ",round(lcm_smry$Conc[lcm_smry$Sample_ID==sample],digits=0)),pos=4)
    }
    
    
    add_annotation(tree=tree,
                   details=details_targ_full,
                   matrices,
                   annot_function=function(tree,details,matrices,node) {
                     add_var_col(tree,
                                 details,
                                 matrices,
                                 node,
                                 var_field = sample,
                                 pval_based=FALSE,
                                 lwd = 5,
                                 colours=colour.scale,
                                 scale_muts_to_branch=scale_muts_to_branch)
                   }
    )
    if(plot_cell_frac) {
      add_annotation_targeted(sample,
                              tree=tree,
                              details=details_targ_full,
                              matrices=matrices,
                              annot_function=function(node,sample,tree,details,matrices,cex=0.6) {
                                node_cell_frac=get_node_cell_frac(node,gsub("_comb","",sample),tree,details,matrices)
                                info=get_edge_info(tree,details,node)
                                if(!is.na(node_cell_frac) & any(post.prob.mat[info$idx.in.details,sample]>prob_threshold_to_include)) {
                                  text(info$x,info$yb,round(node_cell_frac,digits=3),cex = cex,col="black",font=2)
                                }
                              })
    }
    if(plot_donut) {
      #Detect which nodes to plot donuts for - do it for branches with a mean clean.post.prob of >0.5
      nodes_to_check=unique(tree$edge[,2])[!unique(tree$edge[,2])%in%1:length(tree$tip.label)]
      nodes_to_include=sapply(nodes_to_check,function(node) {
        if(mean(post.prob.mat[details_targ$node==node,sample])>0.5){return(node)}else{return(NA)}
      })
      nodes_to_include<-nodes_to_include[!is.na(nodes_to_include)]
      
      if(donut_info=="cell_frac") {
        #Iterate through these nodes and plot the donuts
        for(node in nodes_to_include) {
          print(node)
          data=node_lineage_loss(node=node,sample=sample,tree = tree,details = details_targ,matrices=matrices,boot_straps = 10000,CI=CI,return_ancestral_cell_frac = TRUE)
          #Make the pie chart for plotting on the node
          df2<-data%>%dplyr::select(-median)%>%gather(key="category",value="count")
          df2$ymax = df2$count
          df2$ymin = c(0, head(df2$ymax, n=-1))
          df2=rbind(df2,data.frame(category="lineages_lost",count=(1-df2$ymax[df2$category=="upper_CI"]),ymax=1,ymin=(df2$ymax[df2$category=="upper_CI"])))
          df2$prop=df2$ymax-df2$ymin
          
          #Get the node co-ordinates
          info=get_edge_info(node=node,tree=tree,details=details_targ)
          #Plot the pie chart
          #plotDonut(info$x,mean(c(info$yb,info$yt)),median=data$median,radius=radius,col=c( "#8D8DCB" ,"#C6C6E5", "#FFFFFF"),prop=df2$prop,border="black")
          
          median_only_prop=c(data$median,1-data$median)
          
          plotDonut(info$x,mean(c(info$yb,info$yt)),radius=radius,col=c( "#08306B" ,"#FFFFFF"),prop=median_only_prop,border="black",plotPie = TRUE)
        }
      } else if(donut_info=="lineages_lost") {
        
        #Iterate through these nodes and plot the donuts
        for(node in nodes_to_include) {
          print(node)
          data=node_lineage_loss(node=node,sample=sample,tree = tree,details = details_targ,matrices=matrices,boot_straps = 10000,CI=CI,return_ancestral_cell_frac = FALSE)
          #Make the pie chart for plotting on the node
          df2<-data%>%select(-median)%>%gather(key="category",value="count")
          df2$ymax = df2$count
          df2$ymin = c(0, head(df2$ymax, n=-1))
          df2=rbind(df2,data.frame(category="lineages_lost",count=(1-df2$ymax[df2$category=="upper_CI"]),ymax=1,ymin=(df2$ymax[df2$category=="upper_CI"])))
          df2$prop=df2$ymax-df2$ymin
          
          #Get the node co-ordinates
          info=get_edge_info(node=node,tree=tree,details=details_targ)
          #Plot the pie chart
          plotDonut(info$x,info$yb,median=data$median,radius=radius,col=c( "#8D8DCB" ,"#C6C6E5", "#FFFFFF"),prop=df2$prop,border="black",plotPie = TRUE)
        }
      }
    }
  }
  )
}

#The "squash tree" function to cut the tree at any given node height.  Tree structure is maintained, but edge lengths are shortened.
squash_tree=function(tree,cut_off=50) {
  require(phytools)
  tree$edge.length[nodeHeights(tree)[,1]>=cut_off] <-0 #Any edge that starts at or above the cut-off -> 0
  idxs_to_squash=which(nodeHeights(tree)[,1]<=cut_off & nodeHeights(tree)[,2]>cut_off) #Find the edges that start below the cut-off but end-up above it
  new_edge_lengths=cut_off - nodeHeights(tree)[idxs_to_squash,1] #work-out the edge lengths that these should be such that they finish at the cut-off
  tree$edge.length[idxs_to_squash] <- new_edge_lengths #Assign these edge.lengths to the edges
  return(tree)
}

#Function to calculate vaf, but avoiding dividing by zero by skipping sites with depth of 0.
calculate_vaf=function(NV,NR){
  NR[NR==0]<-1
  vaf=NV/NR
  return(vaf)
}

#Function for getting trees from different stages of development - set the edges to the appropriate length and
#this function prunes the branches that should now be removed
prune_tree_of_zero_tips = function(tree) {
  current_tree=NULL
  for(i in 1:10) {
    print(i)
    if(is.null(current_tree)){current_tree<-tree} else {current_tree<-tree_pruned}
    current_tree$tip.label<-1:length(current_tree$tip.label)
    tips_to_remove=NULL
    for(j in 1:length(current_tree$tip.label)) {
      private_branch_length=current_tree$edge.length[current_tree$edge[,2]==j]
      if(private_branch_length>0) {
        next
      } else {
        ancestor=current_tree$edge[current_tree$edge[,2]==j,1]
        if(ancestor==current_tree$edge[1,1]) {
          ancestral_branch_length<-1
        } else {
          ancestral_branch_length=current_tree$edge.length[current_tree$edge[,2]==ancestor]
        }
        all_daughters=current_tree$edge[current_tree$edge[,1]==ancestor,2]
        daughter_branch_lengths<-current_tree$edge.length[current_tree$edge[,2]%in%all_daughters]
        if(ancestral_branch_length==0|all(daughter_branch_lengths==0)) {
          tips_to_remove=c(tips_to_remove,j)
        }
      }
    }
    if(is.null(tips_to_remove)) {stop(return(current_tree))}
    tree_pruned=drop.tip(current_tree,trim.internal=FALSE,current_tree$tip.label[tips_to_remove])
  }
}

#The "clean up" function to make it fit better with phylogeny
clean_up_post=function(post.prob,details,tree) { #post.prob is the named vector of posterior probs for a given sample (i.e. vector)
  
  #Define a "node map" to allow tracking to true parent nodes in polytomous tree that is structured as dichotomous
  node_map=list()
  for(node in unique(details$node)) {node_map[[node]] = get_direct_daughters(node,tree_targ)}
  
  #Set up the output prob vector to be the same as the input, but can be edited
  clean.post=post.prob
  
  for(mut in details$mut_ref) {
    #Get the node for the current mutation
    allocated_node=details$node[details$mut_ref==mut]
    info=get_edge_info(tree,details,allocated_node)
    #The the daughter node numbers
    daughters=get_direct_daughters(allocated_node,tree=tree)
    #Get the ancestral node number
    ancestor=which(unlist(lapply(node_map,function(x) allocated_node%in%x)))
    
    #Other mutations on the current branch
    other_muts_on_same_branch=details$mut_ref[details$node==allocated_node & details$mut_ref!=mut]
    mean_current_branch=mean(post.prob[other_muts_on_same_branch])
    mean_ancestral_branch_prob=mean(post.prob[details$mut_ref[details$node==ancestor]])
    daughter_probs=post.prob[details$mut_ref[details$node%in%daughters]]
    
    #Remove artefacts, i.e. mutations that don't make sense from rest of phylogeny
    if((length(ancestor)==0|mean_ancestral_branch_prob<0.1) & all(daughter_probs<0.1)) {clean.post[mut] <- 0}
    
    #Boost post.prob of mutations that have confident downstream probabilities
    if((length(ancestor)==0|mean_ancestral_branch_prob>0.5) & (mean_current_branch>0.5|length(other_muts_on_same_branch)==0) & any(daughter_probs >0.8)) {clean.post[mut] <-(post.prob[mut]+1)/2}
  }
  return(clean.post)
}


##Tim's Binomial mixture model functions

## Expectation step
estep = function(x,size,p.vector,prop.vector,ncomp){
  ## p.vector = vector of probabilities for the individual components
  ## prop.vector = vector of proportions for the individual components
  ## ncomp = number of components
  p.mat_estep = matrix(0,ncol=ncomp,nrow=length(x))
  for (i in 1:ncomp){
    p.mat_estep[,i]=prop.vector[i]*dbinom(x,size,prob=p.vector[i])
  }
  norm = rowSums(p.mat_estep) ## normalise the probabilities
  p.mat_estep = p.mat_estep/norm
  LL = sum(log(norm)) ## log-likelihood
  
  ## classification of observations to specific components (too crude?)
  which_clust = rep(1,length(x))
  if(ncomp>1){
    which_clust = apply(p.mat_estep, 1, which.max)
  }
  
  list("posterior"=p.mat_estep,
       "LL"=LL,
       "Which_cluster"=which_clust)
}

## Maximisation step
mstep = function(x,size,e.step){
  # estimate proportions
  prop.vector_temp = colMeans(e.step$posterior)
  # estimate probabilities
  p.vector_temp = colSums(x/size*e.step$posterior) / colSums(e.step$posterior)
  
  list("prop"=prop.vector_temp,
       "p"=p.vector_temp)  
}

## EM algorithm
em.algo = function(x,size,prop.vector_inits,p.vector_inits,maxit=5000,tol=1e-6,nclust){
  ## prop.vector_inits =  initial values for the mixture proportions
  ## p.vector_inits =  initial values for the probabilities
  
  # Initiate EM
  flag = 0
  e.step = estep(x,size,p.vector = p.vector_inits,prop.vector = prop.vector_inits,ncomp=nclust)
  m.step = mstep(x,size,e.step)
  prop_cur = m.step[["prop"]]
  p_cur = m.step[["p"]]
  cur.LL = e.step[["LL"]]
  LL.vector = e.step[["LL"]]
  
  # Iterate between expectation and maximisation steps
  for (i in 2:maxit){
    e.step = estep(x,size,p.vector = p_cur,prop.vector = prop_cur,ncomp=nclust)
    m.step = mstep(x,size,e.step)
    prop_new = m.step[["prop"]]
    p_new = m.step[["p"]]
    
    LL.vector = c(LL.vector,e.step[["LL"]])
    LL.diff = abs((cur.LL - e.step[["LL"]]))
    which_clust = e.step[["Which_cluster"]]
    # Stop iteration if the difference between the current and new log-likelihood is less than a tolerance level
    if(LL.diff < tol){ flag = 1; break}
    
    # Otherwise continue iteration
    prop_cur = prop_new; p_cur = p_new; cur.LL = e.step[["LL"]]
    
  }
  if(!flag) warning("Didn’t converge\n")
  
  BIC = log(length(x))*nclust*2-2*cur.LL
  AIC = 4*nclust-2*cur.LL
  list("LL"=LL.vector,
       "prop"=prop_cur,
       "p"=p_cur,
       "BIC"=BIC,
       "AIC"=AIC,
       "n"=nclust,
       "Which_cluster"=which_clust)
}

binom_mix = function(x,size,nrange=1:3,criterion="BIC",maxit=5000,tol=1e-6){
  ## Perform the EM algorithm for different numbers of components
  ## Select best fit using the Bayesian Information Criterion (BIC)
  ## or the Akaike information criterion (AIC)
  i=1
  results = list()
  BIC_vec = c()
  AIC_vec = c()
  
  for (n in nrange){
    ## Initialise EM algorithm with values from kmeans clustering
    init = kmeans(x/size,n)
    prop_init = init$size/length(x)
    p_init = init$centers
    
    results[[i]] = em.algo(x,size,prop.vector_inits = prop_init,p.vector_inits=p_init,nclust=n,maxit,tol)
    BIC_vec = c(BIC_vec,results[[i]]$BIC)
    AIC_vec = c(AIC_vec,results[[i]]$AIC)
    i=i+1
  }
  if (criterion=="BIC"){
    results[[which.min(BIC_vec)]]$BIC_vec=BIC_vec
    return(results[[which.min(BIC_vec)]])
  }
  if (criterion=="AIC"){
    return(results[[which.min(AIC_vec)]])
  }
}

#This is the function to use for males (as is the case for both foetuses)
calculate_cell_frac=function(NV,NR) {
  #Remove 0's from the depth, to avoid dividing by 0
  NR[NR==0]<-1
  #Get vectors to select out the autosomal and XY chromosomal mutations
  XY_muts=grepl("X",rownames(NR))|grepl("Y",rownames(NR))
  #For autosomal mutations, cell frac = NV/ (NR/2)
  cell_frac=NV/(NR/2)
  #For XY mutations, cell frac = NV/NR -> replace these accordingly
  cell_frac[XY_muts,]<-(NV[XY_muts,])/(NR[XY_muts,])
  #Cell frac cannot be greater than 1, therefore if comes out as > 1 (which can happen when dividing the NR by 2), coerce to 1
  cell_frac[cell_frac>1]<-1
  return(cell_frac)
}
