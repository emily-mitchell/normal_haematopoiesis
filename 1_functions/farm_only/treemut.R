# Note you need to be in the foleder containing this script along with the scripts treemut.c, treemut.Darwin.so and treemut.so in order to source.
require("ape")
require("parallel")

load_lib=function(){
  if(system("uname",intern = TRUE)=="Darwin"){
    dyn.load("treemut.Darwin.so")
  }else{
    dyn.load("treemut.so")
  }
}
load_lib()


reconstruct_genotype_summary=function(phylo
){
  dat=phylo$edge
  samples=phylo$tip.label
  N=length(samples)
  zeros=rep(0,N)
  profile=sapply(1:length(samples),function(i){tmp=zeros;tmp[i]=1;paste(tmp,collapse="")})
  df=data.frame(profile=profile,edge_length=phylo$edge.length[1:N])
  #Create empty list where each element will correspond to a node
  muts=lapply(1:dim(dat)[1],function(i) c())
  map_node=match(1:max(dat[,2]),dat[,2])
  #Here we loop through the tips (samples) and get all the ancestral nodes for each tip.
  #Then for each node in the list we add the tip to that nodes list of tips, ultimately giving us a list of samples that share the node. 
  for(i in 1:N){
    parents=get_ancestral_nodes(i,edge = dat,exclude_root=TRUE)
    for(j in parents){
      muts[[map_node[j]]]=append(muts[[map_node[j]]],i)
    }
  }
  profiles=sapply(muts,function(x){tmp=zeros;tmp[x]=1;paste(tmp,collapse="")})
  df=data.frame(profile=profiles,edge_length=phylo$edge.length,stringsAsFactors=FALSE)
  df=add_derived_profile_info(df,phylo$tip.label)
  list(df=df,samples=phylo$tip.label)
}

get_ancestral_nodes=function(node,edge,exclude_root=TRUE){
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

add_derived_profile_info=function(profile_df,samples=sprintf("s%s",0:(nchar(profile_df$profile[1])-1))){
  base=rep(0,nchar(profile_df$profile[1]))
  samples_private=sapply(1:length(base),function(i){this=base;this[i]=1;paste(this,collapse="")})
  missing_samples=setdiff(samples_private,profile_df$profile)
  if(length(missing_samples)>0){
    profile_df=rbind(profile_df,data.frame(profile=missing_samples,edge_length=0))
  }
  
  profile_df$mut_count=nchar(profile_df$profile)-nchar(gsub("1","",profile_df$profile))
  profile_df$profile_int=lapply(profile_df$profile,split_profile)
  profile_df$id=1:dim(profile_df)[1]
  profile_df$label=profile_df$profile
  idx=which(profile_df$mut_count==1)
  profile_df$label[idx]=sapply(idx,function(i) samples[which(profile_df$profile_int[[i]]==1)])
  profile_df
}


assign_to_tree=function(mtr,dep,df,error_rate=rep(0.01,dim(mtr)[2]),maxits=5){
  if(!is.matrix(mtr)){
    nm=names(mtr)
    mtr=matrix(mtr,nrow=1)
    colnames(mtr)=nm
    dep=matrix(dep,nrow=1)
    colnames(dep)=nm
  }
  p.err=error_rate#c(error_rate,1e-10)
  p.err=p.err[match(df$samples,colnames(mtr))]
  mtr=mtr[,df$samples]
  dep=dep[,df$samples]
  if(!is.matrix(mtr)){
    nm=names(mtr)
    mtr=matrix(mtr,nrow=1)
    colnames(mtr)=nm
    dep=matrix(dep,nrow=1)
    colnames(dep)=nm
  }
  tree_genotypes=do.call("rbind",df$df$profile_int)
  el=rep(1,length(df$df$edge_length))
  loglik=rep(NA,maxits)
  n=dim(mtr)[1]
  for(i in 1:maxits){
    ol=el
    lik=get_likelihood_mtr_C(mtr,dep,tree_genotypes,el,p.error = p.err)
    edge_ml=apply(lik,1,which.max)
    #n=length(edge_ml)
    p=exp(lik-lik[(edge_ml-1)*n+1:n])##Substract max to just control range of lik... comes out in the wash later.
    p=p/rowSums(p)
    loglik[i]=sum(p*lik)
    el=colSums(p,na.rm=T)
    epsilon=sum(abs(el-ol))/dim(mtr)[1]
    cat("delta edge length=",epsilon,"\n")
    cat("Loglik=",loglik[i],"\n")
    if(epsilon<0.01/dim(mtr)[1]){
      break
    }
  }
  
  df$df$expected_edge_length=el
  df$df$edge_length=sapply(1:length(df$df$edge_length),function(i){length(which(edge_ml==i))})
  
  p_else_where=1-p[(edge_ml-1)*n+1:n]
  cat("Finished assigning mutations\ncalculating pvalues\n")
  pval=rep(NA,length(edge_ml))
  for(i in 1:length(pval)){
    pval[i]= get_mutation_assignment_pval(df$df,edge_ml[i],mtr[i,],dep[i,],p.err)
    if(i %% 1000 == 0){
      cat("On",i," of ",length(pval),"\n")
    }
  }
  
  list(df=df,lik=lik,summary=data.frame(edge_ml=edge_ml,pval=pval,p_else_where=p_else_where),p=p)
}

assign_to_tree_bin=function(reported_genotype,fpr=0.05,fnr=0.05,df,maxits=5){
  if(!is.matrix(reported_genotype)){
    nm=names(reported_genotype)
    reported_genotype=matrix(reported_genotype,nrow=1)
    colnames(reported_genotype)=nm
  }
  reported_genotype=reported_genotype[,df$samples]
  if(!is.matrix(reported_genotype)){
    nm=names(reported_genotype)
    reported_genotype=matrix(reported_genotype,nrow=1)
    colnames(reported_genotype)=nm
  }
  
  tree_genotypes=do.call("rbind",df$df$profile_int)
  el=rep(1,length(df$df$edge_length))
  loglik=rep(NA,maxits)
  n=dim(reported_genotype)[1]
  for(i in 1:maxits){
    ol=el
    lik=get_likelihood_reported_genotype(reported_genotype,tree_genotypes,el,fpr,fnr)
    edge_ml=apply(lik,1,which.max)
    #n=length(edge_ml)
    p=exp(lik-lik[(edge_ml-1)*n+1:n])##Substract max to just control range of lik... comes out in the wash later.
    p=p/rowSums(p)
    loglik[i]=sum(p*lik)
    el=colSums(p,na.rm=T)
    epsilon=sum(abs(el-ol))/dim(mtr)[1]
    cat("delta edge length=",epsilon,"\n")
    cat("Loglik=",loglik[i],"\n")
    if(epsilon<0.01/dim(mtr)[1]){
      break
    }
  }
  
  df$df$expected_edge_length=el
  df$df$edge_length=sapply(1:length(df$df$edge_length),function(i){length(which(edge_ml==i))})
  
  p_else_where=1-p[(edge_ml-1)*n+1:n]
  cat("Finished assigning mutations\ncalculating pvalues\n")
  pval=rep(NA,length(edge_ml))
  for(i in 1:length(pval)){
    pval[i]= 1 ##get_mutation_assignment_pval(df$df,edge_ml[i],mtr[i,],dep[i,],p.err)
    if(i %% 1000 == 0){
      cat("On",i," of ",length(pval),"\n")
    }
  }
  
  list(df=df,lik=lik,summary=data.frame(edge_ml=edge_ml,pval=pval,p_else_where=p_else_where),p=p)
}

get_likelihood_mtr_C=function(mtr,depth,geno,el,p.error=rep(0.01,dim(mtr)[2])){
  if(dim(mtr)[2]!=dim(geno)[2]){
    stop("error: dimension mismatch betwee tree genotypes and mtr")
  }
  nmuts=dim(mtr)[1]
  nsamp=dim(mtr)[2]
  nbranch=dim(geno)[1]
  res = .C("likelihood_mtr", 
           mtr=as.integer(mtr),
           depth=as.integer(depth),
           geno=as.integer(geno),
           el=as.double(el),
           p.error=as.double(p.error),
           nmuts=as.integer(nmuts),
           nsamp=as.integer(nsamp),
           nbranch=as.integer(nbranch),
           lik=double(nmuts*nbranch)
  )
  matrix(res[["lik"]],ncol=nbranch)
}

##STUB to be implemented....
get_likelihood_reported_genotype=function(reported_genotype,tree_genotypes,el,fpr,fnr){
  lik=apply(reported_genotype,1,function(rg){
    sapply(1:dim(tree_genotypes)[1],function(x){
      
    })
    
    
  })
  
  
}


get_mutation_assignment_pval=function(df,i,mtr,dep,p.error){
  profile_int=df$profile_int[[i]]
  idx.inside=which(profile_int==1)
  idx.outside=which(profile_int==0)
  V=0.5-p.error
  
  ##Calculate probability for observing <= total MTR at variant sites give VAF=0.5.  This is complicated by allowing different error rates for each sample (otherwise just pbinom)
  if(length(idx.inside)>0){
    ##We breakdown errors into unique categories to make calculation more tractable - could further discretise.
    v=unique(V[idx.inside])
    if(length(v)>4){
      stop("Too many error categories in get_pval: not tractable!")
    }
    mtri=mtr[idx.inside]
    depi=dep[idx.inside]
    vi=V[idx.inside]
    mtrii=sapply(v,function(v) sum(mtri[vi==v]))
    depii=sapply(v,function(v) sum(depi[vi==v]))
    p1=cpv(sum(mtrii),depii,v,lower.tail = TRUE)
    d1=1
  }else{
    p1=1
    d1=0
  }
  ##Calculate probability for observing >= total MTR at non-variant sites give VAF=0.  This is complicated by allowing different error rates for each sample (otherwise just pbinom)
  if(length(idx.outside)>0){
    v=unique(p.error[idx.outside])
    if(length(v)>4){
      stop("Too many error categories in get_pval: not tractable!")
    }
    mtri=mtr[idx.outside]
    depi=dep[idx.outside]
    vi=p.error[idx.outside]
    mtrii=sapply(v,function(v) sum(mtri[vi==v]))
    depii=sapply(v,function(v) sum(depi[vi==v]))
    p2=cpv(sum(mtrii),depii,v,lower.tail = FALSE)
    d2=1
  }else{
    p2=1
    d2=0
  }
  ##We are interested in whether either of the tests fails so we combine using conservative bonferroni.  
  ##We are performing 2 tests.  same as p.adjust(c(p1,p2),method = "bonferroni")
  pv=min((d1+d2)*min(c(p1,p2)),1)
  pv
}


cpv=function(mtrtot,depth,probs,lower.tail){
  ##Calculates the probability that there are a total of MTR mutant reads across N bins each with depth d_i and mutant prob p_i
  ##Does this by explicitly calculating  P(M | {d_i},{p_i)=sum_m P(mtr_i=m | d_i, p_i)*P(M-m | )
  n=length(depth)
  idx=order(depth)
  probs=probs[idx]
  depth=depth[idx]
  flag=ifelse(lower.tail,1,0)
  resk=list(mtrtot,probs,depth,n,lower.tail)
  res=.C("cumulate_binomial_pval",as.integer(round(mtrtot)),as.integer(round(depth)),as.double(probs),as.integer(n),as.integer(lower.tail),p=double(1))
  res[["p"]]
}

split_profile=function(profile){
  as.integer(unlist(strsplit(profile,split="")))
}

##Generate a random tree
generate_random_tree=function(nsamples){
  tree=rtree(nsamples)
  tree=bind.tree(tree,read.tree(text="(zeros:0);"))
  tree$edge.length=ceiling(200*tree$edge.length) ##Scales the #Muts per edge to the scale 0-200
  tree
}

#get_tree_summary_df=function(tree){
#

#####Crude Simulation code
simulate_reads_from_tree=function(df,avgdepth,n_artifacts=0,p.error=0.01){
  if(length(p.error)>1){
    stop("p.error must be scalar here")
  }
  nsamples=length(df$samples)-1
  #df=reconstruct_genotype_summary(tree)
  N=length(df$df$edge_length)
  dat=lapply(1:N,function(i) get_simulated_reads(df$df$profile_int[[i]],df$df$edge_length[i],avgdepth,p.error=p.error))
  #browser()
  
  ml=do.call("c",lapply(1:N,function(i) rep(i,df$df$edge_length[i])))
  if(n_artifacts>0){
    ##Add in some artefacts that aren't simulated according to the tree topology.
    NART=n_artifacts
    dat2=lapply(1:NART,function(i) get_simulated_reads(ifelse(runif(nsamples+1)<0.3,1,0),1,avgdepth,p.error=p.error))
    dat=c(dat,dat2)
    ml=c(ml,rep(NA,NART))
  }
  mtr=do.call("rbind",lapply(dat,function(x) x$mtr))
  depth=do.call("rbind",lapply(dat,function(x) x$depth))
  colnames(mtr)=df$samples
  colnames(depth)=df$samples
  p.error=rep(p.error,length(df$samples))
  idx.zeros=match("zeros",df$samples)
  mtr[,idx.zeros]=0
  depth[,idx.zeros]=10
  p.error[idx.zeros]=1e-6 ##This will force the inference code to keep the zeros branch as a zero length outgroup. 
  list(tree=tree,mtr=mtr,depth=depth,edge=ml,df=df,p.error=p.error)
}

get_simulated_reads=function(geno,n,depth,p.error=0.01){
  if(n==0){
    emptymat=matrix(1,nrow=0,ncol=length(geno))
    return(list(mtr=emptymat,depth=emptymat))
  }
  geno=rep(geno,n)
  depth=rpois(length(geno),depth)
  mtr=rbinom(length(depth),depth,prob=ifelse(geno==0,p.error,0.5-p.error))
  list(mtr=matrix(mtr,nrow=n,byrow = TRUE),
       depth=matrix(depth,nrow=n,byrow = TRUE)
  )
}

run_sim=function(nsamp,depth){
  dat=simulate_reads_from_tree(nsamp,depth)
  res=assign_to_tree(dat$mtr,dat$depth,dat$df,error_rate=dat$p.error)##
  list(edge_length_orig=dat$df$df$edge_length,
       edge_length_inferred=res$df$df$edge_length,
       expected_edge_length_inferred=res$df$df$expected_edge_length,
       edge_idx_orig=dat$edge,
       edge_idx_ml=res$summary$edge_ml)
}

plot_sim_result=function(res,dat){
  sim=list(edge_length_orig=dat$df$df$edge_length,
           edge_length_inferred=res$df$df$edge_length,
           expected_edge_length_inferred=res$df$df$expected_edge_length,
           edge_idx_orig=dat$edge,
           edge_idx_ml=res$summary$edge_ml)
  
  par(mfrow=c(2,2))
  mismatch_prop=length(which(sim$edge_idx_orig!=sim$edge_idx_ml))/length(sim$edge_idx_orig)
  plot(sim$edge_length_orig,sim$edge_length_inferred,xlab="Original Edge Length",ylab="Inferred Edge Length",
       main=sprintf("Hard Assigned edge length vs Orig (sd=%3.2f)\n Mismatch Proportion=%5.4f",
                    sd(sim$edge_length_orig-sim$edge_length_inferred),mismatch_prop)
  )
  abline(a=0,b=1)
  plot(sim$edge_length_orig,sim$expected_edge_length_inferred,xlab="Original Edge Length",ylab="Expected Edge Length",
       main=sprintf("Expected edge length (Soft Assigned) \nvs Orig (sd=%3.2f)",sd(sim$edge_length_orig-sim$expected_edge_length_inferred))
  )
  abline(a=0,b=1)
  simdf=data.frame(deviation=sim$expected_edge_length_inferred-sim$edge_length_orig,
                   edge_length_orig=sim$edge_length_orig)
  simdf=simdf[order(simdf$edge_length_orig),]
  
  with(simdf,plot(edge_length_orig,deviation,
       xlab="Original Edge Length",ylab="Deviation",
       main="Deviation vs Original Edge Length"))
  ##Add LOESS to see if there is any systematic patterns
  loess_fit <- with(simdf,loess(deviation ~ edge_length_orig))
  pl=predict(loess_fit,se = T)
 
  plu=pl$fit+1.96*pl$se
  pll=pl$fit-1.96*pl$se
  
  
  shade_between(simdf$edge_length_orig,plu,pll,adjustcolor("blue",0.3))
  with(simdf,lines(edge_length_orig, pl$fit, col = "blue",lwd=2))
  
  plot(simdf$edge_length_orig,abs(simdf$deviation),
       xlab="Original Edge Length",ylab="Absolute Deviation",
       main="Absolute Deviation vs Original Edge Length")
  
  
  
}

shade_between=function(x,y1,y2,color){
  polygon(c(x, rev(x), x[1]), c(y1, rev(y2), y1[1]), 
          col = color,border = color) 
}


get_newick2=function(profile_df){
  root=1:nchar(profile_df$profile[1])
  children=get_children2(profile_df,root)
  idx.root=which(profile_df$profile==paste(rep(1,length(root)),collapse=""))
  if(length(idx.root)==0){
    mut=0
  }else{
    mut=profile_df$edge_length[idx.root]
  }
  newick_string=sprintf("(%s)R:%s;",
                        paste(
                          sapply(children,
                                 function(child){
                                   s1=get_newick_string2(child,profile_df)
                                   mut=profile_df$edge_length[child]
                                   sprintf("%s:%s",s1,mut)
                                 }
                          ),
                          collapse=","),mut)
  newick_string
}
get_children2=function(profile_df,idx,DEPTH_LIMIT=100){
  k=length(idx)
  children=c()
  #i=1
  depth_check=0
  #browser()
  #cat("testing")
  while(length(idx)>0 && depth_check<DEPTH_LIMIT){
    #cat(idx,"\n")
    depth_check=depth_check+1
    #Find sub-profiles
    outside=sapply(1:dim(profile_df)[1],function(j) sum(profile_df$profile_int[[j]][-idx]))
    inside=sapply(1:dim(profile_df)[1], function(j) sum(profile_df$profile_int[[j]][idx]))
    idx_tree_violations=which(outside>0 & inside>0 & inside!=length(idx))
    if(length(idx_tree_violations)>0){
      
      cat(idx,"\n")
      cat("violated at ",idx_tree_violations,"\n")
      stop("ERROR:Data is non-tree like!")
    }
    descendents=which(outside==0 & inside<k)##Second clause stops parent being included
    if(length(descendents)==0){
      return(children)
    }
    ##Now identify direct descendents
    idx_largest_sub_clade=descendents[which.max(profile_df$mut_count[descendents])]
    idx.child=which(profile_df$profile_int[[idx_largest_sub_clade]]==1)
    children=append(children,profile_df$id[idx_largest_sub_clade])
    idx=setdiff(idx,idx.child)
  }
  return(children)
}
get_newick_string2=function(node,profile_df){
  idx=which(profile_df$profile_int[[node]]==1)
  mut_count=profile_df$edge_length[node]
  children=get_children2(profile_df,idx)
  if(length(children)>0){
    sprintf("(%s)%s",
            paste(
              sapply(children,
                     function(child){
                       s1=get_newick_string2(child,profile_df)
                       mut=profile_df$edge_length[child]
                       sprintf("%s:%s",s1,mut)
                     }
              ),
              collapse=","),profile_df$label[node])
  }else{
    sprintf("%s",profile_df$label[node],mut_count)
  }
}


split_profile=function(profile){
  as.integer(unlist(strsplit(profile,split="")))
}


get_tree_from_df=function(df){
  read.tree(text=get_newick2(df$df))
}

