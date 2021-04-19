library(devtools)
library(ape)
library(phytools)
library(MCMCglmm)
library(phangorn)
library(spam)
library(INLA)
library(phylodyn)
library(ggtree)

#Source PC's functions
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
      if (curr.dist == 0) {curr.dist <- 0.01} # So that no edge lengths are zero
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

# Throw away function to compare methods
generate.diagnostic.plots <- function(tree) {
  #Make the tree dichotomous - required for phylodyn
  tree<- multi2di(tree)
  
  # Generate different versions of ultrametric tree
  tree.nnls <- force.ultrametric(tree, method="nnls")
  tree.extend <- force.ultrametric(tree, method="extend")
  tree.bespoke <- make.ultrametric.tree(tree)
  
  par(mfrow=c(2,2))
  plot(tree, show.tip.label = FALSE, main="Unadjusted")
  plot(tree.nnls, show.tip.label = FALSE, main="NNLS")
  plot(tree.extend, show.tip.label = FALSE, main="Extend")
  plot(tree.bespoke, show.tip.label = FALSE, main="Bespoke")
  
  # Generate population size trajectories
  tree.BNPR <- BNPR(tree)
  tree.BNPR.nnls <- BNPR(tree.nnls)
  tree.BNPR.extend <- BNPR(tree.extend)
  tree.BNPR.bespoke <- BNPR(tree.bespoke)
  
  plot_BNPR(tree.BNPR, col = "red", main="Population trajectory: unadjusted", xlab = "Time", ylab = "Relative population size")
  plot_BNPR(tree.BNPR.nnls, col = "red", main="Population trajectory: NNLS method", xlab = "Time", ylab = "Relative population size")
  plot_BNPR(tree.BNPR.extend, col = "red", main="Population trajectory: extend method", xlab = "Time", ylab = "Relative population size")
  plot_BNPR(tree.BNPR.bespoke, col = "red", main="Population trajectory: custom method", xlab = "Time", ylab = "Relative population size")
  
}
