#Tree functions written by Nick Williams / Mike Spencer Chapman

#get_ancestral_nodes
#get_idx_for_node
#get_edge_info
#add_annotation
#add_binary_proportion
#add_simple_labels
#add_vaf
#plotBars
#plotPie
#plot_tree_vaf
#plot_tree_labels_genes
#get_color_pch_df
#get_qdf
#plot_tree_labels


require("RColorBrewer")

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

get_idx_for_node=function(details,node){
  which(details$node==node)
}

get_edge_info=function(tree,details,node){
  y=get_y_range(tree,node)
  x=get_x_range(tree,node)
  idx=get_idx_for_node(details,node)
  samples=get_samples_in_clade(node,tree)
  list(yb=y[1],yt=y[2],x=x[1],xm=x[2],idx.in.details=idx,samples=samples)
}

add_annotation=function(tree,details,matrices,annot_function){
  N=dim(tree$edge)[1]
  lapply(1:N,function(i) annot_function(tree,details,matrices,tree$edge[i,2]))
}

add_binary_proportion=function(tree,##<< enhanced phylo returned from plot_tree
                               details,##<< dataframe with summary details of mutations mapped to tree together with EDGE_IDX - index of mutation in edges matrix
                               matrices,##<< a list of matrices parallel to details with columns named by tip labels
                               node,
                               bfield,
                               b.add.line=TRUE,
                               b.add.text=FALSE,
                               ...
                               ){
  ##Get all the detail about the edge coords + idx in detail
  info=get_edge_info(tree,details,node)
   
  bdat=details[[bfield]][info$idx]
  if(is.null(bdat) || class(bdat)!="logical"){
    stop("Error in provided bfield (does it exist and is it boolean?)")
  }
  pass=sum(bdat,na.rm=TRUE)
  fail=sum(!bdat,na.rm=TRUE)
  tot=pass+fail
  ycross=info$yb+(fail/tot)*(info$yt-info$yb)
  ##Could add in a third category NA
  #missing=sum(is.na(bdat))
  if(b.add.line){
    arrows(y0=info$yb,y1=ycross,x0=info$x,length = 0,col="black",lend=1,...)
    arrows(y0=ycross,y1=info$yt,x0=info$x,length = 0,col="red",lend=1,...)
  }
  if(b.add.text){
    text(y=ycross,x=info$x,label=pass,pos = 4,offset = 0,...)
  }
}

add_simple_labels=function(tree,##<< enhanced phylo returned from plot_tree
                           details,##<< dataframe with summary details of mutations mapped to tree using the "node" column -> tree$edge[,2]
                           matrices,##<< a list of matrices parallel to details with columns named by tip labels
                           node,##<< Node (see details)
                           query.field,##<< Name of column in details to query against
                           query.allowed.df,##<< Values of query field which should be annotated. data.frame value,col,pch columns.
                           label.field,##<< Name of column in details specifying the label text.
                           cex.label=1,
                           b.add.label=TRUE,
                           b.add.marker=TRUE,
                           ... ##<< paremeters for points (not color)
){
  info=get_edge_info(tree,details,node)
  idx=info$idx[which(details[[query.field]][info$idx] %in% query.allowed.df$value)]
  if(length(idx)>50){
    stop("too many variants to annotate")
  }
  query.value=details[[query.field]][idx]
  idx.match=match(query.value,query.allowed.df$value)
  cols=query.allowed.df$col[idx.match]
  pch=query.allowed.df$pch[idx.match]
  
  vlabels=details[[label.field]][idx]
  ## spread out
  N=length(idx)
  ##Vertical offset so that labels sit slightly above the markers.
  voffset=0.0075*(par("usr")[4]-par("usr")[3])
  if(N>0){
    yd=info$yt-info$yb
    if(N==1){
      y=0.5*(info$yb+info$yt)
    }else{
      y=seq(info$yb+(1/(N+1))*yd,info$yt-(1/(N+1))*yd,length.out = N)
    }
    if(b.add.marker){
      points(rep(info$x,N),y,col=cols,pch=pch,...)
    }
    if(b.add.label){
      text(rep(info$x,N),y+voffset,labels = vlabels,pos = 2,offset = 0,cex=cex.label)
    }
  }
  list(node=node,value=query.value)
}

add_vaf=function(tree,##<< enhanced phylo returned from plot_tree
                 details,##<< dataframe with summary details of mutations mapped to tree together with EDGE_IDX - index of mutation in edges matrix
                 matrices,##<< a list of matrices parallel to details with columns named by tip labels
                 node,
                 samples=NULL,
                 b.plot.bars=TRUE,
                 lwd.rect=1,
                 min.depth=1,
                 vc.field,
                 vc.df,
                 filter.on=NULL,
                 ...
){
  ##Get all the detail about the edge coords + idx in detail
  info=get_edge_info(tree,details,node)
  #browser()
  ##Can do additional filter e.g. missense 
  #if(!is.null(filter.on)){
  #info$idx=info$idx,3)#info$idx.in.details[which(details$VC=="missense")]
  # }
  
  if(length(info$idx)==0){
    return(NULL)
  }
  if(b.plot.bars){
    plotF=plotBars
  }else{
    plotF=plotPie
  }
  if(is.null(samples)){
    samples=info$samples
  }
  cat(length(info$idx),"\n")
  if(length(samples)>1){
    if(length(info$idx)>1){
      df=data.frame(mtr=rowSums(matrices$mtr[info$idx,samples],na.rm = TRUE),
                    dep=rowSums(matrices$dep[info$idx,samples],na.rm = TRUE),stringsAsFactors = FALSE)
    }else{
      df=data.frame(mtr=sum(matrices$mtr[info$idx,samples],na.rm = TRUE),
                    dep=sum(matrices$dep[info$idx,samples],na.rm = TRUE),stringsAsFactors = FALSE)
    }
  }else{
    df=data.frame(mtr=matrices$mtr[info$idx,samples],
                  dep=matrices$dep[info$idx,samples],stringsAsFactors = FALSE)
  }
  df=cbind(df,details[info$idx,])
  df=df[which(df$dep>=min.depth),]
  df$vaf=df$mtr/df$dep
  df=df[which(!is.na(df$vaf)),]
  N=dim(df)[1]
  if(N==0){
    return(df)
  }
  
  
  df=df[order(df$vaf),]
  yd=info$yt-info$yb
  
  
  if(N==1){
    y=0.5*(info$yb+info$yt)
    width=yd
  }else{
    y1=seq(info$yb,info$yt,length.out = N+2)
    #Don't use the ends..
    y=y1[2:(N+1)]
    width=y[2]-y[1]
  }
  
  if(!b.plot.bars){
    r=0.8  ##r>0.5 will cause likely overlap problems
  }else{
    r=0.4
  }
  #arrows(x0=info$x-w,x1=info$x-w+df$vaf*2*w,y0=y,lend=1,length=0,col="black")
  #arrows(x0=info$x-w+df$vaf*2*w,x1=info$x+w,y0=y,lend=1,length=0,col="grey")
  for(i in 1:N){
    vaf=min(df$vaf[i],0.999)
    if(is.na(vaf)){
      plotF(x=info$x,y = y[i],radius=r,col=c("lightgray","lightgray"),prop=c(0,1),border="lightgray",width=width)
    }else{
      plotF(x=info$x,y = y[i],radius = r,col=c("black","white"),prop = c(vaf,1-vaf),width = width)
    }
  }
  if( !b.plot.bars){
    
    return(df)
  }
  ##Now check to see if we need to highlight
  ##Test if VAF is significantly > 0.05 or significantly < 0.45
  ##Can also do a binomial test...
  MTR=sum(df$mtr)
  DEP=sum(df$dep)
  min.mean.vaf=0.45
  z=binom.test(MTR,DEP,alternative = "less",p=min.mean.vaf)
  z2=binom.test(MTR,DEP,alternative = "greater",p=0.05)
  z$p.value=max(z$p.value,z2$p.value)
  txt=gsub("^0\\.",".",sprintf("%3.2f",MTR/DEP))
  if(z$p.value<0.05){
    if(z$p.value<0.05/dim(tree$edge)[1]){
      border.color="red"
    }else{
      border.color="blue"
    }
  }else{
    border.color="darkgrey"
  }
  
  rect(xleft=info$x-r,xright=info$x+r,ybottom=y[1]-width/2,ytop=y[N]+width/2,border=border.color,lwd=lwd.rect)
  if(border.color!="darkgrey"){
    text(txt,x=info$x,y=y[1]+0.3*(y[N]-y[1]),col="black",cex=0.6)
  }
  arrows(x0=info$x,y0=info$yb,y1=info$yt,lwd=0.5,col="black",length=0,lend=2)
  
  df
  
}


plotBars=function(x,y,radius,col,prop,border="black",width=1){
  #cat(prop,"\n")
  if(width<2){
    arrows(x0 = x-radius,y0=y,x1=x-radius+2*radius*prop[1],col="darkgrey",lend=2,length=0)
    arrows(x0 = x-radius+2*radius*prop[1],y0=y,x1=x-radius+2*radius,col=rgb(0.98,0.98,0.98),lend=2,length=0)
  }else{
    rect(xleft = x-radius,xright =x-radius+2*radius*prop[1],ybottom = y-width/2,ytop=y+width/2,border = NA,col="darkgrey")
    rect(xleft =  x-radius+2*radius*prop[1],xright =x+radius,ybottom = y-width/2,ytop=y+width/2,border = NA,col=rgb(0.98,0.98,0.98))
  }
  1
}

plotPie=function(x,y,radius,col,prop,llwd=0.5,border="black",width=NA){
  lims=par("usr")
  as=dev.size()
  asr=as[1]/as[2]
  yscale=asr*(lims[4]-lims[3])/(lims[2]-lims[1])
  prop=prop/sum(prop)
  cutpoint=c(0,cumsum(prop)*2*pi)
  
  N=2*pi/0.05
  n=ceiling(N*diff(cutpoint)/(2*pi))
  d=diff(cutpoint)/n
  if(length(prop)>1){
    for(i in 2:length(cutpoint)){
      polygon(x+c(radius*cos(seq(cutpoint[i-1],cutpoint[i],d[i-1])),0),
              y+yscale*c(radius*sin(seq(cutpoint[i-1],cutpoint[i],d[i-1])),0),
              border=border,col=col[i-1],lwd=llwd)
    }
  }else{
    i=2
    polygon(x+c(radius*cos(seq(cutpoint[i-1],cutpoint[i],d[i-1]))),
            y+yscale*c(radius*sin(seq(cutpoint[i-1],cutpoint[i],d[i-1]))),
            border=border,col=col[1],lwd=llwd)
  }
  yscale
}

plot_tree_vaf=function(tree,details,matrices,samples=NULL){
  tree=plot_tree(tree)
  res=add_annotation(tree,details,matrices,
                     function(tree,details,matrices,node){
                       add_vaf(tree,details,matrices,node,samples=samples,b.plot.bars = FALSE)
                       }
                     )
  ##Post process res to add legend
  
}


plot_tree_labels_genes=function(tree,details,query.field="GENE",label.field="GENE",genes=c("JAK2","CBL","TET2","DNMT3A"),cex.label=1){
  qdf=data.frame(value=genes,col=rainbow(length(genes)),pch=19)
  plot_tree_labels(tree,details,
                   query.allowed.df = qdf,
                   query.field=query.field,
                   label.field=label.field,
                   cex.label=cex.label)
}



##Gets unique colour pch combos and returns in dataframe with columns "col" and "pch"
get_color_pch_df=function(n){
  pch.list=c(18,17,16,15,0:6)
  if(n>length(pch.list)*8){
    stop("Too many colours requested")
  }
  cols=rep(RColorBrewer::brewer.pal(8,"Set1"),times=length(pch.list))
  pch=rep(pch.list,each=8)
  data.frame(col=cols,pch=pch,stringsAsFactors = FALSE)[1:n,]
  
}

get_qdf=function(values){
  if(length(values)>length(unique(values))){
    stop("get_qdf: please provide values without duplication")
  }
  cbind(data.frame(value=values,stringsAsFactors = FALSE),
        get_color_pch_df(length(values)))
}

plot_tree_labels_consequence=function(tree,details,consequences,
                                      query.allowed.df=get_qdf(consequences),
                          query.field="VC",
                          label.field="GENE",
                          cex.label=1){
  ##qdf=get_qdf(consequences)
  plot_tree_labels(tree,details,
                   query.allowed.df = query.allowed.df,
                   query.field=query.field,
                   label.field=label.field,
                   cex.label=cex.label)
}


plot_tree_labels=function(tree,details,
                          query.field="VC",
                          query.allowed.df=data.frame(value=c("nonsense","frameshift"),
                                                      col=c("red","black"),pch=c(17,18)
                                                      ),
                          label.field="GENE",
                          cex.label=1){
  res=add_annotation(tree,
                     details,list(),
                     function(tree,details,matrices,node){
                       add_simple_labels(tree,details,matrices,node,
                                         query.field =query.field,
                                         query.allowed.df = query.allowed.df,
                                         label.field = label.field,
                                         cex.label =cex.label)})
  with(query.allowed.df,legend("topleft",legend=value,col=col,pch=pch))
}


plot_tree_vaf=function(tree,details,matrices,samples=NULL,b.plot.bars =TRUE,filter.on=NULL){
  res=add_annotation(tree,details,matrices,
                     function(tree,details,matrices,node){
                       add_vaf(tree,details,matrices,node,samples=samples,b.plot.bars = b.plot.bars,filter.on=filter.on)})
}

#add_vaf_bar
#add_vof_pie
#add_label





