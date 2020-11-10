#' @importFrom Matrix t
#' @importFrom SummarizedExperiment assay
.collate_se_matrices <- function(x, assays) {
    targets1 <- lapply(assays, FUN=assay, x=x)
    lapply(targets1, t)
}

#' @importFrom SingleCellExperiment reducedDim altExp
.collate_sce_matrices <- function(x, assays=NULL, dimreds=NULL, altexps=NULL, altexp.assay="logcounts") {
    targets2 <- targets3 <- list()

    if (!is.null(dimreds)) {
        targets2 <- lapply(dimreds, FUN=reducedDim, x=x)
    }

    if (!is.null(altexps)) {
        targets3 <- lapply(altexps, FUN=altExp, x=x)
        altexp.assay <- rep(altexp.assay, length.out=length(targets3))
        targets3 <- mapply(FUN=assay, x=targets3, i=altexp.assay, SIMPLIFY=FALSE)
        targets3 <- lapply(targets3, t)
    }

    c(targets2, targets3)
}
