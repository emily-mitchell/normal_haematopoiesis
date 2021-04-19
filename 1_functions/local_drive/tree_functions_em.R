#Tree functions adapted by me and functions for making ultrametric tree

#plot_tree_em
#find.distance
#length.normalise
#make.ultrametric.tree
#generate.bespoke.plots

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