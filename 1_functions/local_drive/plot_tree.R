#Tree functions written by Nick Williams

#set_cedge
#get_height
#elbow
#elbowv
#set_tree_coords
#plot_tree (I have used an amended version in functions_for_tree_visualisation_em)
#get_all_node_children
#get_node_children
#get_samples_in_clade
#get_y_range
#get_x_range



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

