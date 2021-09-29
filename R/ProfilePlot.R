#' Plot the phylogenetic profiles
#'
#' Plot the phylogenetic profiles along with the phylogenetic tree.
#' This function is based on `heatmap.2()` from `gplots` package.
#' @import gplots
#' @param profile a data frame or matrix contains the phylogenetic profiles.
#' Each column is a profile. The row names must match the tip labels of the tree.
#' @param ... some other parameters available for `heatmap.2()`, such as `main`, `xlab`, `labRow`...
#' @param dend the dendrogram that will be plotted on the left side. See `TreeToDen()`.
#' @import gplots
#' @examples
#' d = TreeToDend(rtree(50))
#' p = matrix(sample(0:1,100, replace=T), nrow=50)
#' rownames(p) = labels(d)
#' ProfilePlot(p, d, main="Plot of two profiles")
#' @export

ProfilePlot = function(profile,dend, ... ){
    if(!is.data.frame(profile) & !is.matrix(profile)){
        stop("the profile should be a data frame or matrix")
    }
    if (!setequal(rownames(profile), labels(dend))){
        stop("the rownames of profile are not matched with the tree")
        }
    mat = as.matrix(profile[labels(dend),])
    gplots::heatmap.2(mat,Rowv=dend,dendrogram="row",col =c("white","black"),
              sepcolor=c("white"),sepwidth=c(0.05,0.05),key=FALSE,Colv = F,trace = "none",tracecol = "yellow",colsep = 1:ncol(mat),rowsep=1:nrow(mat),
              cexRow=0.5,cexCol = 0.8, ...)

}
