#' Intersect two graphs
#'
#' Intersect two graphs by taking the product of their edge weights.
#'  
#' @param ... Any number of \link{graph} objects for the same set and order of nodes.
#' @param mixing Numeric vector of length equal to the number of graphs, specifying the mixing weights for the edge weights.
#' @param nominal Numeric scalar specifying the scaling factor to compute the nominal weight.
#'
#' @return A \link{graph} object containing the intersection of \code{g1} and \code{g2}.
#'
#' @details
#' The idea of taking an intersection of a weighted graph is based on the intersection of simplicial sets in the UMAP algorithm.
#' For each edge that exists in either graph, we compute the product of the weights across all graphs and assign that value as the edge weight in the output graph.
#' This means that edges in the output only have high weight if they are present and highly weighted in all graphs - hence, an intersection.
#'
#' This approach would make for a very sparse graph if the product was taken directly.
#' To maintain some connectivity, edges that exist in one graph but not the other are assigned nominal weights in the latter to ensure that the product is not zero.
#' The nominal weight for each graph is defined as the product of its smallest non-zero edge weight and \code{nominal}.
#' Decreasing this value will yield a more conservative intersection and a less connected graph, usually manifesting as smaller clusters after application of community detection algorithms.
#' 
#' By default, \code{mixing} is a vector of length equal to the numbeer of graphs, containing values of 1.
#' This means that edge weights from each graph in \code{...} contribute equally to the product.
#' However, it is possible to increase the contribution of some of the graphs by supplying a higher \code{mixing} values for those graphs.
#'
#' Unweighted graphs are supported and are considered to have edge weights of 1.
#' 
#' @author Aaron Lun
#'
#' @examples
#' library(scran)
#'
#' mat1 <- matrix(rnorm(10000), ncol=20)
#' chosen <- 1:250
#' mat1[chosen,1] <- mat1[chosen,1] + 10
#' g1 <- buildSNNGraph(mat1, d=NA, transposed=TRUE)
#' clusters1 <- igraph::cluster_walktrap(g1)$membership
#' table(clusters1, chosen=mat1[,1] > 5)
#'
#' # Pretending we have some other data for the same cells, e.g., ADT.
#' mat2 <- matrix(rnorm(10000), ncol=20)
#' chosen <- c(1:125, 251:375)
#' mat2[chosen,2] <- mat2[chosen,2] + 10
#' g2 <- buildSNNGraph(mat2, d=NA, transposed=TRUE)
#' clusters2 <- igraph::cluster_walktrap(g2)$membership
#' table(clusters2, mat2[,2] > 5)
#' 
#' # Intersecting the graphs and clustering:
#' gcom <- intersectGraphs(g1, g2)
#' clustersC <- igraph::cluster_walktrap(gcom)$membership
#' table(clustersC, mat1[,1] > 5)
#' table(clustersC, mat2[,2] > 5)
#' @export
#' @importFrom igraph E ends make_graph E<-
#' @importFrom BiocGenerics unique rbind match
#' @importFrom S4Vectors DataFrame
#' @importFrom scuttle .unpackLists
intersectGraphs <- function(..., mixing=NULL, nominal=1e-6) {
    g.list <- list(...)

    p.list <- w.list <- vector("list", length(g.list))
    for (i in seq_along(g.list)) {
        curg <- g.list[[i]]
        e <- ends(curg, E(curg))
        df <- DataFrame(first=e[,1], second=e[,2])
        w <- E(curg)$weight

        keep <- w > 0 # protect against zero-weight edges.
        p.list[[i]] <- df[keep,,drop=FALSE]
        w.list[[i]] <- w[keep]
    }

    U <- unique(do.call(rbind, p.list))
    for (i in seq_along(p.list)) {
        w <- w.list[[i]]
        if (is.null(w)) {
            w <- rep(1, nrow(p.list[[i]]))
        }
        minw <- min(w[w > 0]) * nominal

        m <- match(U, p.list[[i]])    
        w <- w[m]
        w[is.na(w)] <- minw
        w.list[[i]] <- w
    }

    if (is.null(mixing)) {
        mixing <- rep(1, length(g.list))
    }
    mixing <- mixing/min(mixing) # based on how UMAP does it.
    w.list2 <- mapply(`^`, w.list, mixing, SIMPLIFY=FALSE)

    output <- make_graph(rbind(U$first, U$second), directed=FALSE)
    E(output)$weight <- Reduce("*", w.list2)
    output
}
