#' Convert tree to dendrogram
#'
#' Convert the non-ultrametric tree to a dendrogram for `ProfilePlot()` to plot the profiles using
#' the `chronos()` function from `ape` package
#'
#' @param tree an object of class "phylo"
#' @return a dendrogram as an input for `ProfilePlot()` function
#' @import ape
#' @examples
#' #library(ape)
#' # create a random 50-tip tree
#' t = rtree(50)
#' plot(t)
#' d = TreeToDend(t)
#' plot(d)
#' @export


TreeToDend = function(phytree){
    ## transfer .nex format tree to a dendrogram
    phytree = ape::chronos(phytree)
    dend=stats::as.dendrogram(ape::as.hclust.phylo(phytree))
    return(dend)
}
