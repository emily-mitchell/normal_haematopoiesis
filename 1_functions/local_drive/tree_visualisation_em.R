# Tree visualisation functions amended by Emily Mitchell

#get_ancestral_nodes 
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