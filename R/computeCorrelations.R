#' Compute correlations between modes
#'
#' Compute Spearman correlations between two sets of features, using data collected for the same cells in different modalities.
#'
#' @inheritParams findTopCorrelations
#'
#' @return
#' A DataFrame where each row corresponds to a pair of features in \code{x} and \code{y}.
#' (If \code{y=NULL}, each pair corresponds to a pair of features in \code{x}.)
#' This contains the following fields:
#' \itemize{
#' \item \code{feature1}, the name (character) or row index (integer) of each feature in \code{x}.
#' \item \code{feature2}, the name (character) or row index (integer) of one of the top correlated features to \code{feature1}.
#' This is another feature in \code{x} if \code{y=NULL}, otherwise it is a feature in \code{y}.
#' \item \code{rho}, the Spearman rank correlation for the current pair of \code{feature1} and \code{feature2}.
#' \item \code{p.value}, the approximate p-value associated with \code{rho} under the null hypothesis that the correlation is zero.
#' \item \code{FDR}, the adjusted p-value.
#' }
#' The rows are sorted by \code{feature1} and then \code{p.value}.
#'
#' @details
#' If \code{block} is specified, correlations are computed separately for each block of cells.
#' For each feature pair, the reported \code{rho} is set to the average of the correlations across all blocks.
#' If \code{equiweight=FALSE}, the average is weighted by the number of cells in each block.
#'
#' Similarly, the p-value corresponding to each correlation is computed separately for each block and then combined across blocks with Stouffer's method.
#' More specifically, combining is done using the one-sided p-values for both signs of the correlation, and the smaller p-value is taken (and multiplied by 2).
#' This ensures that a low p-value can only be achieved if the blocks agree in the sign.
#' If \code{equiweight=FALSE}, each per-block p-value is weighted by the number of cells.
#'
#' @author Aaron Lun
#' @examples
#' library(scuttle)
#' sce1 <- mockSCE()
#' sce1 <- logNormCounts(sce1)
#'
#' sce2 <- mockSCE(ngenes=10) # pretend this is protein data.
#' sce2 <- logNormCounts(sce2)
#'
#' output <- computeCorrelations(sce1, sce2)
#' output
#'
#' @seealso
#' \code{\link{findTopCorrelations}}, to avoid computing correlations for all pairs of features when \code{y} has many rows.
#'
#' @name computeCorrelations
NULL

######################
##### Internals ######
######################

#' @importFrom Matrix t
#' @importFrom beachmat rowBlockApply
#' @importFrom S4Vectors DataFrame
#' @importFrom stats p.adjust
#' @importFrom metapod parallelStouffer
#' @importFrom BiocParallel SerialParam bpstart bpstop
#' @importFrom scuttle .bpNotSharedOrUp .subset2index
#' @importFrom scran rhoToPValue
.compute_all_correlations <- function(x, y, subset.cols=NULL, block=NULL, equiweight=TRUE, use.names=TRUE, BPPARAM=SerialParam()) {
    if (!.bpNotSharedOrUp(BPPARAM)) {
        bpstart(BPPARAM)
        on.exit(bpstop(BPPARAM))
    }

    if (!is.null(subset.cols)) {
        subset.cols <- .subset2index(subset.cols, x, byrow=FALSE)
        x <- x[,subset.cols,drop=FALSE]
        y <- y[,subset.cols,drop=FALSE]
        block <- block[subset.cols]
    }

    if (is.null(block)) {
        by.block <- NULL
        nblocks <- ncol(x)
    } else {
        by.block <- split(seq_len(ncol(x)), block)
        nblocks <- lengths(by.block)
    }

    if (is.null(y)) {
        alt <- x
    } else {
        alt <- y
    }
    rho <- rowBlockApply(x, FUN=.compute_exact_rho, other=alt, by.block=by.block, grid=TRUE, BPPARAM=BPPARAM)
    rho <- do.call(mapply, c(list(FUN=rbind, SIMPLIFY=FALSE), rho)) # still fragmented from the blockApply.
    mean.rho <- .compute_mean_rho(rho, nblocks, equiweight)

    if (is.null(y)) {
        names1 <- names2 <- rownames(x)
        nother <- nrow(x)
    } else {
        names1 <- rownames(x)
        names2 <- rownames(y)
        nother <- nrow(y)
    }

    self <- rep(seq_len(nrow(x)), nother)
    other <- rep(seq_len(nother), each=nrow(x))
    df <- DataFrame(feature1=self, feature2=other, rho=as.vector(mean.rho))
    df <- .fill_names(df, use.names, names1, names2)

    p.values <- mapply(FUN=rhoToPValue, rho=rho, n=nblocks, SIMPLIFY=FALSE)

    if (length(nblocks)==1L) {
        up.value <- p.values[[1]]$positive
        down.value <- p.values[[1]]$negative
    } else {
        # Combining all the one-sided p-values together.
        if (equiweight) {
            weights <- NULL
        } else {
            weights <- nblocks
        }

        up.values <- lapply(p.values, "[[", i="positive")
        up.value <- parallelStouffer(up.values, weights=weights)$p.value
        down.values <- lapply(p.values, "[[", i="negative")
        down.value <- parallelStouffer(down.values, weights=weights)$p.value
    }

    df$p.value <- pmin(up.value, down.value, 0.5)*2
    if (is.null(y)) {
        df <- df[df$feature1!=df$feature2,,drop=FALSE]
    }
    df$FDR <- p.adjust(df$p.value, method="BH")

    o <- order(df$feature1, df$p.value)
    df[o,,drop=FALSE]
}

#' @importFrom stats cor
#' @importFrom beachmat rowBlockApply
.compute_exact_rho <- function(x, other, by.block) {
    FUN <- function(block2, primary) {
        cor(
            as.matrix(t(primary)), 
            as.matrix(t(block2)), 
            method="spearman"
        )
    }

    if (is.null(by.block)) {
        output <- rowBlockApply(other, FUN=FUN, primary=x, grid=TRUE)
        output <- list(do.call(cbind, output))
    } else {
        output <- vector("list", length(by.block))
        for (n in seq_along(by.block)) {
            chosen <- by.block[[n]] 
            subout <- rowBlockApply(other[,chosen,drop=FALSE], FUN=FUN, grid=TRUE, primary=x[,chosen,drop=FALSE])
            output[[n]] <- do.call(cbind, subout)
        }
    }

    output
}

#' @importFrom S4Vectors DataFrame
#' @importFrom metapod averageParallelStats
.compute_mean_rho <- function(rho, nblocks, equiweight) {
    weights <- NULL
    if (!equiweight) weights <- nblocks
    averageParallelStats(rho, weights=weights)
}

.fill_names <- function(df, use.names, names1, names2) {
    if (use.names) {
        if (!is.null(names1)) {
            df$feature1 <- names1[df$feature1]
        }
        if (!is.null(names2)) {
            df$feature2 <- names2[df$feature2]
        }
    }
    df
}

#######################
##### S4 methods ######
#######################

#' @export
#' @rdname computeCorrelations
setGeneric("computeCorrelations", function(x, y, ...) standardGeneric("computeCorrelations"))

#' @export
#' @rdname computeCorrelations
setMethod("computeCorrelations", "ANY", .compute_all_correlations)

#' @export
#' @rdname computeCorrelations
setMethod("computeCorrelations", "SummarizedExperiment", function(x, y, use.names=TRUE, ..., assay.type="logcounts") {
    formatted <- .format_xy_for_correlations(x, y, use.names=use.names, assay.type=assay.type)
    .compute_all_correlations(formatted$x, formatted$y, use.names=formatted$use.names, ...)
})

#' @importFrom SummarizedExperiment rowData assay
.format_xy_for_correlations <- function(x, y, use.names=TRUE, assay.type="logcounts") {
    x0 <- assay(x, assay.type)

    y0 <- y
    if (is(y, "SummarizedExperiment")) {
        y0 <- assay(y, assay.type)
    }

    if (is.character(use.names)) {
        use.names <- rep(use.names, length.out=2)
        if (!is.na(use.names[1])) {
            rownames(x0) <- rowData(x)[,use.names[[1]]]
        }
        if (!is.na(use.names[2]) && is(y, "SummarizedExperiment")) {
            rownames(y0) <- rowData(y)[,use.names[[2]]]
        }
        use.names <- TRUE
    }

    list(x=x0, y=y0, use.names=use.names)
}


