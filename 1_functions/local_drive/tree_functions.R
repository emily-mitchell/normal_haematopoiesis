#Tree functions (some adapted for current work by Emily Mitchell)

#Authors: Nick Williams, Mike Spencer Chapman, Peter Campbell 

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
#set_cedge
#get_height
#elbow
#elbowv
#set_tree_coords
#get_all_node_children
#get_node_children
#get_samples_in_clade
#get_y_range
#get_x_range
#nodeHeights (function from phytools)
#plot_node_number
#highlight_nodes
#boxtext (online function - used in add_simple_labels_line)
#add_simple_labels
#add_simple_labels_line
#add_simple_line
#plot_tree_lines
#add_simple_labels_only
#plot_tree_labels
#plot_tree
#add_heatmap
#plot_tree_em
#find.distance
#length.normalise
#make.ultrametric.tree
#generate.bespoke.plots
#nodeheight
#getAncestors
#getTips
#correct_edge_length
#get_subset_tree
#get_corrected_tree
#get_mut_burden
#get_mut_burden_stats
#get_summary_stats




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

set_cedge=function(parent,tree){
  for(child in get_node_children(parent,tree)){
    #cat("\nsetting parent - child ",parent,child,"\n")
    child.idx=which(tree$edge[,2]==child)
    parent.idx=which(tree$edge[,2]==parent)
    if(length(parent.idx)==0){
      pedge=0
    }else{
      pedge=tree$cedge[parent.idx]
    }
    #cat("before:",length(which(tree$cedge>0)),"\n")
    tree$cedge[child.idx]=pedge+tree$edge.length[child.idx];
    #cat("after:",length(which(tree$cedge>0)),"\n")
    #cat(parent,child,":",tree$cedge)
    tree=set_cedge(child,tree)
  }
  tree
}

##Not very efficient -- recursively calculates height as average of children's height.
get_height=function(tree,node){
  #if(is.null(tree$height)){
  #  tree$height=rep(NA,1:length(tree$edge.length))
  #}
  N=length(tree$tip.label)
  if(node<=N){
    return(node)#tree$height[which(tree$edge[,2]==node)]=node
  }else{
    children=get_node_children(node,tree)
    return(mean(sapply(children,function(child) get_height(tree,child))))
  }
}

set_height=function(tree){
  tree$height_end=sapply(tree$edge[,2],function(i) get_height(tree,i))
  tree$height_start=tree$height[match(tree$edge[,1],tree$edge[,2])]
  N=length(tree$tip.label)
  root=N+1
  idx=which(tree$edge[,1]==root)
  tree$height_start[idx]=mean(tree$height_end[idx])
  tree
}

##horizontal elbow
elbow=function(x0,x1,y0,y1,...){ 
  arrows(x0=x0,y0=y0,y1=y1,length=0,...)
  arrows(x0=x0,x1=x1,y0=y1,length=0,...)
}
##vertical elbow
elbowv=function(x0,x1,y0,y1,...){ 
  #browser()
  arrows(x0=x0,x1=x1,y0=y0,length=0,...)
  arrows(x0=x1,y0=y0,y1=y1,length=0,...)
  
}

set_tree_coords=function(atree){
  ##get the cumulative distance from the root.
  tree=atree
  tree$cedge=rep(0,length(tree$edge.length))
  N=length(tree$tip.label)
  root=N+1
  tree=set_cedge(root,tree)
  tt=set_height(tree)
  atree$coords=data.frame(a0=tt$cedge-tt$edge.length,a1=tt$cedge,
                          b0=tt$height_start,b1=tt$height_end,stringsAsFactors = FALSE)
  if(!is.null(atree$color)){
    atree$coords$color=atree$color
  }
  atree
}

# plot_tree=function(tree,direction="down",cex.label=1,offset=0,b_do_not_plot=FALSE,lwd=1,bars=NULL,default_edge_color="darkgrey",ymax=NULL,cex.terminal.dots=0){
#   par(mar=c(1, 1, 1, 3) + 0.1)
#   #browser()
#   if(!(direction %in% c("down","across"))){
#     stop("Unsupported direction provided")
#   }
#   N=length(tree$tip.label)
#   if(is.null(tree$coords)){
#     tree=set_tree_coords(tree)
#   }
#   coords=tree$coords
#   
#   if(direction=="across"){
#     xmax=max(coords$a1)*1.05
#     ymax=max(coords$b1)+1
#     offset=offset*xmax
#   }else{
#     if(is.null(ymax)){
#       ymax=max(coords$a1)*1.05
#     }
#     xmax=max(coords$b1)+1
#     offset=offset*ymax
#   }
#   if(b_do_not_plot){
#     return(tree)
#   }
#   if(is.null(bars)){
#     plot(NULL,axes=FALSE,xlim=c(0-(xmax*0.1),xmax),ylim=c(0-(ymax*0.05),ymax),xlab="",ylab="")
#   }else{
#     
#     plot(NULL,axes=FALSE,xlim=c(0-(xmax*0.1),xmax),ylim=c(0-(ymax*0.15),ymax),xlab="",ylab="")
#   }
#   idx.tip=match(1:N,tree$edge[,2])
#   if(direction=="across"){
#     apply(coords,1,function(x) elbow(x[1],x[2],x[3],x[4]))
#     text(tree$tip.label,x =coords$a1[idx.tip]+offset ,y=coords$b1[idx.tip],cex = cex.label,pos = 4)
#   }else{
#     top=max(coords$a1)
#     ##browser()
#     m=dim(coords)[1]
#     if(is.null(coords$color)){
#       col=rep(default_edge_color,m)
#     }else{
#       col=coords$color
#     }
#     sapply(1:m,function(i) {x=as.numeric(coords[i,1:4]);elbowv(x[3],x[4],top-x[1],top-x[2],col=col[i],lwd=lwd)})
#     if(is.null(tree$tip.color)){
#       tipcol="black"
#     }else{
#       tipcol=tree$tip.color
#     }
#     if(cex.label>0){
#       text(tree$tip.label,y =top-(coords$a1[idx.tip]+offset) ,x=coords$b1[idx.tip],cex = cex.label,pos = 1,col=tipcol)
#     }
#     if(cex.terminal.dots>0){
#       points(y =top-(coords$a1[idx.tip]) ,x=coords$b1[idx.tip],col=tipcol,cex=cex.terminal.dots,pch=19)
#     }
#   }
#   tree$direction=direction
#   tree$top=top
#   scales=c(0,1,10,50,100,200,500,1000,2000,5000)
#   #scale=scales[max(which(ymax/2>scales))]
#   scale=scales[max(which(ymax/5>=scales))]
#   #browser()
#   cat("scale=",scale,"\n")
#   axis(side = 4,at=seq(top,-scale,-scale),label=seq(0,top+scale,scale),las=2)
#   #arrows(x0=length(tree$tip.label)+0.5,y0=0,y1=scale,length=0.1,code=3,angle=90)
#   #text(sprintf("%s Muts",scale),x=length(tree$tip.label)-0.5,y=0.5*scale,pos=4,cex=cex.label,offset=0.1)
#   if(!is.null(bars)){
#     maxbar=max(bars)
#     idx=match(names(bars),tree$tip.label)
#     rect(xleft=idx-0.5,xright=idx+0.5,ybottom = -ymax*0.15,ytop=-ymax*0.15+ymax*0.1*bars/maxbar,col = "grey")
#     
#   }
#   tree
# }


get_all_node_children=function(node,tree){
  children=tree$edge[which(tree$edge[,1]==node),2]
  offspring=children
  for(child in children){
    offspring=c(offspring,get_all_node_children(child,tree))
  }
  offspring
}
get_node_children=function(node,tree){
  tree$edge[which(tree$edge[,1]==node),2]
}

get_samples_in_clade=function(node,tree){
  if(node<=length(tree$tip.label)){
    return(tree$tip.label[node])
  }
  tree$tip.label[intersect(get_all_node_children(node,tree),1:length(tree$tip.label))]
}

get_y_range=function(tree,node){
  idx=which(tree$edge[,2]==node)
  if(length(idx)!=1){
    stop("bad node provided")
  }
  as.numeric(tree$top-tree$coords[idx,c("a1","a0")])
}

get_x_range=function(tree,node){
  idx=which(tree$edge[,2]==node)
  if(length(idx)!=1){
    stop("bad node provided")
  }
  as.numeric(tree$coords[idx,c("b1","b0")])
}

nodeHeights = function (tree, ...) 
{
  if (hasArg(root.edge)) 
    root.edge <- list(...)$root.edge
  else root.edge <- FALSE
  if (root.edge) 
    ROOT <- if (!is.null(tree$root.edge)) 
      tree$root.edge
  else 0
  else ROOT <- 0
  nHeight <- function(tree) {
    tree <- reorder(tree)
    edge <- tree$edge
    el <- tree$edge.length
    res <- numeric(max(tree$edge))
    for (i in seq_len(nrow(edge))) res[edge[i, 2]] <- res[edge[i, 
                                                               1]] + el[i]
    res
  }
  nh <- nHeight(tree)
  return(matrix(nh[tree$edge], ncol = 2L) + ROOT)
}

plot_node_number = function(tree,details,matrices,node,cex=0.4) {
  info=get_edge_info(tree,details,node)
  text(info$x,info$yb,node,cex = cex,col="black",font=2)
}

highlight_nodes=function(tree,details,matrices,node,nodes,...) {
  info=get_edge_info(tree,details,node=node)
  if(node %in% nodes){
    arrows(y0=info$yb,y1=info$yt,x0=info$x,x1=info$x,length=0,col="red",lend=1,...)
  }
}

#Useful function from online for drawing background boxes to the labels (then used in the "add_simple_labels_line" function)
boxtext <- function(x, y, labels = NA, col.text = NULL, col.bg = NA, 
                    border.bg = NA, adj = NULL, pos = NULL, offset = 0.5, 
                    padding = c(0.5, 0.5), cex = 1, font = graphics::par('font')){
  
  ## The Character expansion factor to be used:
  theCex <- graphics::par('cex')*cex
  
  ## Is y provided:
  if (missing(y)) y <- x
  
  ## Recycle coords if necessary:    
  if (length(x) != length(y)){
    lx <- length(x)
    ly <- length(y)
    if (lx > ly){
      y <- rep(y, ceiling(lx/ly))[1:lx]           
    } else {
      x <- rep(x, ceiling(ly/lx))[1:ly]
    }       
  }
  
  ## Width and height of text
  textHeight <- graphics::strheight(labels, cex = theCex, font = font)
  textWidth <- graphics::strwidth(labels, cex = theCex, font = font)
  
  ## Width of one character:
  charWidth <- graphics::strwidth("e", cex = theCex, font = font)
  
  ## Is 'adj' of length 1 or 2?
  if (!is.null(adj)){
    if (length(adj == 1)){
      adj <- c(adj[1], 0.5)            
    }        
  } else {
    adj <- c(0.5, 0.5)
  }
  
  ## Is 'pos' specified?
  if (!is.null(pos)){
    if (pos == 1){
      adj <- c(0.5, 1)
      offsetVec <- c(0, -offset*charWidth)
    } else if (pos == 2){
      adj <- c(1, 0.5)
      offsetVec <- c(-offset*charWidth, 0)
    } else if (pos == 3){
      adj <- c(0.5, 0)
      offsetVec <- c(0, offset*charWidth)
    } else if (pos == 4){
      adj <- c(0, 0.5)
      offsetVec <- c(offset*charWidth, 0)
    } else {
      stop('Invalid argument pos')
    }       
  } else {
    offsetVec <- c(0, 0)
  }
  
  ## Padding for boxes:
  if (length(padding) == 1){
    padding <- c(padding[1], padding[1])
  }
  
  ## Midpoints for text:
  xMid <- x + (-adj[1] + 1/2)*textWidth + offsetVec[1]
  yMid <- y + (-adj[2] + 1/2)*textHeight + offsetVec[2]
  
  ## Draw rectangles:
  rectWidth <- textWidth + 2*padding[1]*charWidth
  rectHeight <- textHeight + 2*padding[2]*charWidth    
  graphics::rect(xleft = xMid - rectWidth/2, 
                 ybottom = yMid - rectHeight/2, 
                 xright = xMid + rectWidth/2, 
                 ytop = yMid + rectHeight/2,
                 col = col.bg, border = border.bg)
  
  ## Place the text:
  graphics::text(xMid, yMid, labels, col = col.text, cex = theCex, font = font, 
                 adj = c(0.5, 0.5))    
  
  ## Return value:
  if (length(xMid) == 1){
    invisible(c(xMid - rectWidth/2, xMid + rectWidth/2, yMid - rectHeight/2,
                yMid + rectHeight/2))
  } else {
    invisible(cbind(xMid - rectWidth/2, xMid + rectWidth/2, yMid - rectHeight/2,
                    yMid + rectHeight/2))
  }    
}

########### Original
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
      #text(rep(info$x,N),y+voffset,labels = vlabels,pos = 2,offset = 0,cex=cex.label)
      boxtext(info$x-1,info$yb+0.5*(info$yt-info$yb),col.bg="white",border.bg="black",padding = c(0.5, 0.5),labels = vlabels,pos=2,cex=cex.label)
    }
  }
  list(node=node,value=query.value)
}

########## Original

add_simple_labels_line=function(tree,##<< enhanced phylo returned from plot_tree
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
  if(length(idx)>1){
    stop("Some branches have multiple variants")
  }
  query.value=details[[query.field]][idx]
  idx.match=match(query.value,query.allowed.df$value)
  cols=query.allowed.df$col[idx.match]
  
  vlabels=details[[label.field]][idx]
  ## spread out
  N=length(idx)
  ##Vertical offset so that labels sit slightly above the markers.
  if(N>0){
    arrows(y0=info$yb,y1=info$yt,x0=info$x,x1=info$x,length=0,col="red",lend=1,lwd=2,lty=1,...)
    if(b.add.label){
      #text(rep(info$x,N),y=info$yb+0.5*(info$yt-info$yb),labels = vlabels,pos = 2,offset = 0.25,cex=cex.label)
      boxtext(info$x-1,info$yb+0.5*(info$yt-info$yb),col.bg="white",border.bg="black",padding = c(0.5, 0.5),labels = vlabels,pos=2,cex=cex.label)
    }
  }
  list(node=node,value=query.value)
}

########################

## To first add lines
add_simple_line=function(tree,##<< enhanced phylo returned from plot_tree
                         details,##<< dataframe with summary details of mutations mapped to tree using the "node" column -> tree$edge[,2]
                         matrices,##<< a list of matrices parallel to details with columns named by tip labels
                         node,##<< Node (see details)
                         query.field,##<< Name of column in details to query against
                         query.allowed.df,##<< Values of query field which should be annotated. data.frame value,col,pch columns.
                         label.field,##<< Name of column in details specifying the label text.
                         cex.label=2,
                         b.add.label=FALSE,
                         b.add.marker=TRUE,
                         col = col,
                         lwd = lwd,
                         lty = lty,
                         ... ##<< paremeters for points (not color)
){
  info=get_edge_info(tree,details,node)
  idx=info$idx[which(details[[query.field]][info$idx] %in% query.allowed.df$value)]
  if(length(idx)>1){
    stop("Some branches have multiple variants")
  }
  query.value=details[[query.field]][idx]
  idx.match=match(query.value,query.allowed.df$value)
  cols=query.allowed.df$col[idx.match]
  
  vlabels=details[[label.field]][idx]
  ## spread out
  N=length(idx)
  ##Vertical offset so that labels sit slightly above the markers.
  if(N>0){
    arrows(y0=info$yb,y1=info$yt,x0=info$x,x1=info$x,length=0,col=col,lend=1,lwd=lwd,lty=lty,...)
    #if(b.add.label){
    #text(rep(info$x,N),y=info$yb+0.5*(info$yt-info$yb),labels = vlabels,pos = 2,offset = 0.25,cex=cex.label)
    #  boxtext(info$x-1,info$yb+0.5*(info$yt-info$yb),col.bg="white",border.bg="white",padding = c(0.5, 0.5),labels = vlabels,pos=2,cex=cex.label)
    # }
  }
  list(node=node,value=query.value)
}


plot_tree_lines=function(tree,details,
                         query.field="VC",
                         type="label",
                         query.allowed.df=data.frame(value=c("nonsense","frameshift")),
                         col = "red",
                         lty=1,
                         lwd=8,
                         label.field="GENE",
                         cex.label=1){
  
  if(type=="label") {res=add_annotation(tree,
                                        details,list(),
                                        function(tree,details,matrices,node){
                                          add_simple_line(tree,details,matrices,node,
                                                          query.field =query.field,
                                                          query.allowed.df = query.allowed.df,
                                                          label.field = label.field,
                                                          cex.label =cex.label)})}
  if(type=="line") {res=add_annotation(tree,
                                       details,list(),
                                       function(tree,details,matrices,node){
                                         add_simple_line(tree,details,matrices,node,
                                                         query.field =query.field,
                                                         query.allowed.df = query.allowed.df,
                                                         label.field = label.field,
                                                         cex.label =cex.label,
                                                         col = col,
                                                         lty=lty,
                                                         lwd=lwd)})}
  # with(query.allowed.df,legend("topleft",legend=value,col=col))
}

add_annotation=function(tree,details,matrices,annot_function,...){
  N=dim(tree$edge)[1]
  lapply(1:N,function(i) annot_function(tree,details,matrices,tree$edge[i,2],...))
}


## To add labels second

add_simple_labels_only=function(tree,##<< enhanced phylo returned from plot_tree
                                details,##<< dataframe with summary details of mutations mapped to tree using the "node" column -> tree$edge[,2]
                                matrices,##<< a list of matrices parallel to details with columns named by tip labels
                                node,##<< Node (see details)
                                query.field,##<< Name of column in details to query against
                                query.allowed.df,##<< Values of query field which should be annotated. data.frame value,col,pch columns.
                                label.field,##<< Name of column in details specifying the label text.
                                cex.label=4,
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
      #boxtext(info$x-1,info$yb+0.5*(info$yt-info$yb),col.bg="white",border.bg="white",padding = c(0.5, 0.5),labels = vlabels,pos=2,cex=cex.label)
      #boxtext(rep(info$x,N),y+voffset,col.bg="",border.bg="",padding = c(0.3, 0.3),labels = vlabels,pos=2,cex=cex.label)
    }
  }
  list(node=node,value=query.value)
}


plot_tree_labels=function(tree,details,
                          query.field="final",
                          type="label",
                          query.allowed.df=data.frame(value= "Final", col = "0", pch = 17),
                          label.field="variant_ID",
                          cex.label=4){
  if(type=="label") {res=add_annotation(tree,
                                        details,list(),
                                        function(tree,details,matrices,node){
                                          add_simple_labels_only(tree,details,matrices,node,
                                                                 query.field =query.field,
                                                                 query.allowed.df = query.allowed.df,
                                                                 label.field = label.field,
                                                                 cex.label =cex.label,
                                          )})}
  
  #with(query.allowed.df,legend("topleft",legend=value,col=col))
}



#### Adapted version of plot_tree

plot_tree=function(tree,direction="down",cex.label=5,offset=0,b_do_not_plot=FALSE,lwd=1.5,bars=NULL,default_edge_color="black",ymax=NULL,cex.terminal.dots=0,vspace.reserve=0){
  par(mar=c(1, 1, 1, 20) + 0.1)
  #browser()
  if(!(direction %in% c("down","across"))){
    stop("Unsupported direction provided")
  }
  N=length(tree$tip.label)
  if(is.null(tree$coords)){
    tree=set_tree_coords(tree)
  }
  coords=tree$coords
  
  if(direction=="across"){
    xmax=max(coords$a1)*1.05
    ymax=max(coords$b1)+1
    offset=offset*xmax
  }else{
    if(is.null(ymax)){
      ymax=max(coords$a1)*1.05
    }
    xmax=max(coords$b1)+1
    offset=offset*ymax
  }
  if(b_do_not_plot){
    return(tree)
  }
  if(is.null(bars)){
    ymin=0-ymax*0.05-vspace.reserve*ymax
    plot(NULL,axes=FALSE,xlim=c(0-(xmax*0.1),xmax),ylim=c(ymin,ymax),xlab="",ylab="")
  }else{
    plot(NULL,axes=FALSE,xlim=c(0-(xmax*0.1),xmax),ylim=c(0-(ymax*0.15),ymax),xlab="",ylab="")
  }
  idx.tip=match(1:N,tree$edge[,2])
  if(direction=="across"){
    apply(coords,1,function(x) elbow(x[1],x[2],x[3],x[4]))
    text(tree$tip.label,x =coords$a1[idx.tip]+offset ,y=coords$b1[idx.tip],cex = cex.label,pos = 4)
  }else{
    top=max(coords$a1)
    ##browser()
    m=dim(coords)[1]
    if(is.null(coords$color)){
      col=rep(default_edge_color,m)
    }else{
      col=coords$color
    }
    sapply(1:m,function(i) {x=as.numeric(coords[i,1:4]);elbowv(x[3],x[4],top-x[1],top-x[2],col=col[i],lwd=lwd)})
    if(is.null(tree$tip.color)){
      tipcol="black"
    }else{
      tipcol=tree$tip.color
    }
    #if(cex.label>0){
    #text(tree$tip.label,y =top-(coords$a1[idx.tip]+offset) ,x=coords$b1[idx.tip],cex = cex.label,pos = 1,col=tipcol)
    # }
    if(cex.terminal.dots>0){
      points(y =top-(coords$a1[idx.tip]) ,x=coords$b1[idx.tip],col=c("darkgrey", "blueviolet","deeppink")[Y_loss], cex=cex.terminal.dots,pch=15)
    }
  }
  tree$direction=direction
  tree$top=top
  #scale =10
  scales=c(0,10,100,200,500,1000,2000,5000)
  scale=scales[max(which(ymax/2>scales))]
  #scale=scales[max(which(ymax/2>=scales))]
  #browser()
  cat("scale=",scale,"\n")
  axis(side = 4,at=seq(top,-scale,-scale),label=seq(0,top+scale,scale),las=2, cex.axis = 5, lwd = 5, lwd.ticks = 5, col = "black")
  #arrows(x0=length(tree$tip.label)+0.5,y0=0,y1=scale,length=0.1,code=3,angle=90)
  #text(sprintf("Mutation number",scale),x=length(tree$tip.label)-0.5,y=0.5*scale,pos=4,cex=cex.label,offset=0.1)
  if(!is.null(bars)){
    maxbar=max(bars)
    idx=match(names(bars),tree$tip.label)
    rect(xleft=idx-0.5,xright=idx+0.5,ybottom = -ymax*0.15,ytop=-ymax*0.15+ymax*0.1*bars/maxbar,col = "black")
    
  }
  tree$ymax=ymax
  tree$vspace.reserve=vspace.reserve
  tree
}

## Add heatmap with additional information under tree
add_heatmap=function(tree,heatmap,heatvals=NULL,border="white",cex.label=2){
  ymax=tree$ymax
  idx=match(colnames(heatmap),tree$tip.label)
  top=-0.05*ymax
  gap=tree$vspace.reserve/dim(heatmap)[1]
  labels=rownames(heatmap)
  for(i in 1:dim(heatmap)[1]){
    bot=top-0.05*ymax
    rect(xleft=idx-0.5,xright=idx+0.5,ybottom = bot,ytop=top,col = heatmap[i,],border=border)
    if(!is.null(heatvals)){
      text(xx=idx,y=0.5*(top+bot),labels = sprintf("%3.2f",heatvals[i,]))
    }
    if(!is.null(labels)){
      text(labels[i],x=-0.5,y=0.5*(top+bot),pos = 2,cex = cex.label)
    }
    top=bot
  }
  tree
}

plot_tree_em = function(tree,direction="down",cex.label=1,offset=0,b_do_not_plot=FALSE,lwd=1,bars=NULL,default_edge_color="black",ymax=NULL,cex.terminal.dots=0){
  par(mar=c(1, 1, 1, 3) + 0.1)
  #browser()
  if(!(direction %in% c("down","across"))){
    stop("Unsupported direction provided")
  }
  N=length(tree$tip.label)
  if(is.null(tree$coords)){
    tree=set_tree_coords(tree)
  }
  coords=tree$coords
  
  if(direction=="across"){
    xmax=max(coords$a1)*1.05
    ymax=max(coords$b1)+1
    offset=offset*xmax
  }else{
    if(is.null(ymax)){
      ymax=max(coords$a1)*1.05
    }
    xmax=max(coords$b1)+1
    offset=offset*ymax
  }
  if(b_do_not_plot){
    return(tree)
  }
  if(is.null(bars)){
    plot(NULL,axes=FALSE,xlim=c(0-(xmax*0.1),xmax),ylim=c(0-(ymax*0.05),ymax),xlab="",ylab="")
  }else{
    
    plot(NULL,axes=FALSE,xlim=c(0-(xmax*0.1),xmax),ylim=c(0-(ymax*0.15),ymax),xlab="",ylab="")
  }
  idx.tip=match(1:N,tree$edge[,2])
  if(direction=="across"){
    apply(coords,1,function(x) elbow(x[1],x[2],x[3],x[4]))
    text(tree$tip.label,x =coords$a1[idx.tip]+offset ,y=coords$b1[idx.tip],cex = cex.label,pos = 4)
  }else{
    top=max(coords$a1)
    ##browser()
    m=dim(coords)[1]
    if(is.null(coords$color)){
      col=rep(default_edge_color,m)
    }else{
      col=coords$color
    }
    sapply(1:m,function(i) {x=as.numeric(coords[i,1:4]);elbowv(x[3],x[4],top-x[1],top-x[2],col=col[i],lwd=lwd)})
    if(is.null(tree$tip.color)){
      tipcol="black"
    }else{
      tipcol=tree$tip.color
    }
    if(cex.label>0){
      text(tree$tip.label,y =top-(coords$a1[idx.tip]+offset) ,x=coords$b1[idx.tip],cex = cex.label,pos = 1,col=tipcol)
    }
    if(cex.terminal.dots>0){
      points(y =top-(coords$a1[idx.tip]) ,x=coords$b1[idx.tip],col=tipcol,cex=cex.terminal.dots,pch=19)
    }
  }
  tree$direction=direction
  tree$top=top
  scales=c(0,1,10,50,100,200,500,1000,2000,5000)
  #scale=scales[max(which(ymax/2>scales))]
  scale=scales[max(which(ymax/5>=scales))]
  #browser()
  cat("scale=",scale,"\n")
  axis(side = 4,at=seq(top,-scale,-scale),label=seq(0,top+scale,scale),las=2)
  #arrows(x0=length(tree$tip.label)+0.5,y0=0,y1=scale,length=0.1,code=3,angle=90)
  #text(sprintf("%s Muts",scale),x=length(tree$tip.label)-0.5,y=0.5*scale,pos=4,cex=cex.label,offset=0.1)
  if(!is.null(bars)){
    maxbar=max(bars)
    idx=match(names(bars),tree$tip.label)
    rect(xleft=idx-0.5,xright=idx+0.5,ybottom = -ymax*0.15,ytop=-ymax*0.15+ymax*0.1*bars/maxbar,col = "grey")
    
  }
  tree
}


# Function to convert a tree to an ultrametric tree

find.distance <- function(tree, from, to) {
  path <- nodepath(tree, from, to)
  res <- 0
  for (j in 2:length(path)) {
    index <- which(tree$edge[,2] == path[j] & tree$edge[,1] == path[j-1], arr.ind = TRUE)
    res <- res + tree$edge.length[index]
  }
  return(res)
}

length.normalise <- function(orig.tree, new.tree, curr.node, remaining.stick) {
  curr.node.children <- unlist(Descendants(orig.tree, curr.node, "children"))
  
  for (j in curr.node.children) {
    index <- which(orig.tree$edge[,2] == j & orig.tree$edge[,1] == curr.node, arr.ind = TRUE)
    
    if (j %in% orig.tree$tip.label) {
      new.tree$edge.length[index] <- remaining.stick
    } else {
      curr.node.tips <- unlist(Descendants(orig.tree, j, "tips"))
      curr.dist <- find.distance(orig.tree, curr.node, j)
      if (curr.dist == 0) {curr.dist <- 0.5} # So that no edge lengths are zero
      desc.lengths <- sapply(curr.node.tips, FUN = find.distance, tree = orig.tree, from = curr.node)
      new.tree$edge.length[index] <- remaining.stick * curr.dist / mean(desc.lengths)
      shorter.stick <- remaining.stick - new.tree$edge.length[index]
      
      # Call as recursive function
      new.tree <- length.normalise(orig.tree, new.tree, j, shorter.stick)
    }
  }
  return(new.tree)
} 

make.ultrametric.tree <- function(tree) {
  root.number <- length(tree$tip.label) + 1
  ultra.tree <- length.normalise(tree, tree, root.number, 1)
  return(ultra.tree)
}

generate.bespoke.plots <- function(tree) {
  # Generate different versions of ultrametric tree
  tree.bespoke <- make.ultrametric.tree(tree)
  
  plot(tree.bespoke, show.tip.label = FALSE)
  
  # Generate population size trajectories
  tree.BNPR.bespoke <- BNPR(tree.bespoke)
  plot_BNPR(tree.BNPR.bespoke)
  
}

nodeheight=function (tree, node, ...) 
{
  if (hasArg(root.edge)) 
    root.edge <- list(...)$root.edge
  else root.edge <- FALSE
  if (root.edge) 
    ROOT <- if (!is.null(tree$root.edge)) 
      tree$root.edge
  else 0
  else ROOT <- 0
  if (!inherits(tree, "phylo")) 
    stop("tree should be an object of class \"phylo\".")
  if (node == (Ntip(tree) + 1)) 
    h <- 0
  else {
    a <- setdiff(c(getAncestors(tree, node), node), Ntip(tree) + 
                   1)
    h <- sum(tree$edge.length[sapply(a, function(x, e) which(e == 
                                                               x), e = tree$edge[, 2])])
  }
  h + ROOT
}

getAncestors=function (tree, node, type = c("all", "parent")) 
{
  if (!inherits(tree, "phylo")) 
    stop("tree should be an object of class \"phylo\".")
  type <- type[1]
  if (type == "all") {
    aa <- vector()
    rt <- Ntip(tree) + 1
    currnode <- node
    while (currnode != rt) {
      currnode <- getAncestors(tree, currnode, "parent")
      aa <- c(aa, currnode)
    }
    return(aa)
  }
  else if (type == "parent") {
    aa <- tree$edge[which(tree$edge[, 2] == node), 1]
    return(aa)
  }
  else stop("do not recognize type")
}

getTips = function(tree,node) {
  require(ape)
  if(node <= length(tree$tip.label)) {
    daughters <- tree$tip.label[node]
  } else {
    daughters <- extract.clade(tree, node = node)$tip.label
  }
  return(daughters)
}

#Note - only apply the "get_edge_from_tree" option if starting from an SNV only tree. Function will assume that all the existing edge length is SNVs.
correct_edge_length = function(node, tree, details, sensitivity_df, include_indels = TRUE, include_SNVs = TRUE, get_edge_from_tree=FALSE) {
  daughters <- getTips(tree = tree, node = node)
  #correct SNVs on edge, or set to 0 if want an indel only tree
  if(include_SNVs == TRUE) {
    if(get_edge_from_tree) {
      nSNV=tree$edge.length[tree$edge[,2]==node]
    } else {
      nSNV = sum(details$node == node & details$Mut_type == "SNV")
    }   		
    all_sens_SNVs <- sensitivity_df[sensitivity_df$Sample %in% daughters,"SNV_sensitivity"]
    branch_SNV_sens = 1 - prod(1-all_sens_SNVs)  
    new_nSNV = nSNV/branch_SNV_sens
  } else {
    new_nSNV <- 0
  }
  #correct INDELs on edge, or set to 0 if want an SNV only tree
  if(include_indels == TRUE) {
    nINDEL = sum(details$node == node & details$Mut_type == "INDEL")
    all_sens_INDELs <- sensitivity_df[sensitivity_df$Sample %in% daughters,"INDEL_sensitivity"]
    branch_INDEL_sens = 1 - prod(1-all_sens_INDELs)
    new_nINDEL = nINDEL/branch_INDEL_sens
  } else {
    new_nINDEL <- 0
  }
  new_edge_length = new_nSNV + new_nINDEL
  return(new_edge_length)
}

get_subset_tree = function(tree, details, v.field = "Mut_type", value = "SNV") {
  get_new_edge_length = function(node, tree, details,v.field,value) {
    sum(details$node == node & details[v.field] == value)
  }
  tree_subset = tree
  tree_subset$edge.length = sapply(tree$edge[,2], get_new_edge_length, tree = tree, details = details,v.field = v.field,value=value)
  return(tree_subset)
}

get_corrected_tree = function(tree, details, sensitivity_df, include_indels = TRUE, include_SNVs = TRUE,get_edge_from_tree=FALSE) {
  tree_c = tree
  tree_c$edge.length = sapply(tree$edge[,2], correct_edge_length, tree = tree, details = details, sensitivity_df = sensitivity_df, include_indels = include_indels, include_SNVs=include_SNVs,get_edge_from_tree=get_edge_from_tree)
  return(tree_c)
}

get_mut_burden = function(tree) {
  mut_burden = nodeHeights(tree)[tree$edge[,2] %in% 1:length(tree$tip.label),2]
}

get_mut_burden_stats = function(tree) {
  mut_burden = get_mut_burden(tree)
  cat(paste("Mean mutation burden is", round(mean(mut_burden),digits = 1),"\n"))
  cat(paste("Range of mutation burden is", round(range(mut_burden)[1],digits = 1),"to",round(range(mut_burden)[2],digits = 1),"\n"))
  cat(paste("Standard deviation of mutation burden is", round(sd(mut_burden),digits = 1),"\n"))
}


get_summary_stats=function(tree,time_points) {
  nodeheights <- nodeHeights(tree)
  sumstats_each_timepoint=lapply(time_points, function(time_point) {
    clades_post_cutoff=tree$edge[,2][nodeheights[,1] < time_point & !nodeheights[,2] < time_point]
    clade_sizes=sapply(clades_post_cutoff,function(node) {length(getTips(tree,node))})
    n_singletons=sum(clade_sizes==1)
    largest_clade=max(clade_sizes)
    return(data.frame(time_point=time_point,n_singletons=n_singletons,largest_clade=largest_clade))
  })
  df=Reduce(cbind,sumstats_each_timepoint)
  return(df)
}

