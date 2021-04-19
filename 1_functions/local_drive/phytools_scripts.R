# nodeHeights function from phytools only

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
