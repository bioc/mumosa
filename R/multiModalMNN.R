#' Multi-modal batch correction with MNNs
#'
#' Perform MNN correction on multi-modal data, based on a generalization of \code{\link{fastMNN}} to multiple feature sets.
#'
#' @param ... One or more \linkS4class{SingleCellExperiment} objects,
#' containing a shared set of alternative Experiments corresponding to different data modalities.
#' Alternatively, one or more lists of such objects.
#' @param batch Factor specifying the batch to which each cell belongs, when \code{...} contains only one SingleCellExperiment object.
#' Otherwise, each object in \code{...} is assumed to contain cells from a single batch.
#' @param which Character vector containing the names of the alternative Experiments to use for correction.
#' Defaults to the names of all alternative Experiments that are present in every object of \code{...}.
#' @param rescale.k Integer scalar specifying the number of neighbors to use in \code{\link{rescaleByNeighbors}}.
#' @param pca.common.args Further arguments to pass to \code{\link{multiBatchPCA}} to control the PCA for all modalities.
#' @param pca.main.args Further arguments to pass to \code{\link{multiBatchPCA}} to control the PCA for the main Experiments.
#' Overrides any arguments of the same name in \code{pca.common.args}.
#' @param pca.alt.args Further arguments to pass to \code{\link{multiBatchPCA}} to control the PCA for each alternative Experiment specified in \code{which}.
#' This should be a list where each entry is named after any alternative Experiment and contains an internal list of named arguments;
#' these override any settings in \code{pca.common.args} in the PCA for the corresponding modality.
#' @param mnn.args Further arguments to pass to \code{\link{reducedMNN}}, controlling the MNN correction.
#' @param BNPARAM A \linkS4class{BiocNeighborParam} object specifying how the nearest neighbor searches should be performed.
#' @param BPPARAM A \linkS4class{BiocParallelParam} object specifying how parallelization should be performed.
#' 
#' @details
#' This function implements a multi-modal MNN correction for SingleCellExperiment inputs where each main and alternative Experiment corresponds to one modality.
#' We perform a PCA within each modality with \code{\link{multiBatchPCA}},
#' rescale the PCs to be of a comparable scale with \code{\link{rescaleByNeighbors}},
#' and finally correct in low-dimensional space with \code{\link{reducedMNN}}.
#' Corrected expression values for each modality are then recovered in the same manner as described for \code{\link{fastMNN}}.
#'
#' Modality-specific arguments can be passed to the PCA step via the \code{pca.common.args}, \code{pca.main.args} and \code{pca.alt.args} arguments.
#' These mirror the corresponding arguments in \code{\link{applyMultiSCE}} - see the documentation for that function for more details.
#' Additional arguments for the MNN step can be passed via \code{mnn.args}.
#' Note that \code{batch} is used across all steps and must be specified as its own argument in the \code{multiModalMNN} function signature.
#'
#' @return
#' A \linkS4class{SingleCellExperiment} of the same structure as that returned by \code{\link{fastMNN}},
#' i.e., with a \code{corrected} entry of corrected low-dimensional coordinates and a \code{reconstructed} assay of corrected expression values.
#' In addition, the \code{\link{altExps}} entries contain corrected values for each data modality used in the correction.
#' 
#' @author Aaron Lun
#' @seealso
#' \code{\link{fastMNN}}, for MNN correction within a single modality.
#'
#' \code{\link{multiBatchPCA}}, to perform a batch-aware PCA within each modality.
#'
#' \code{\link{applyMultiSCE}}, which inspired this interface for Experiment-specific arguments.
#' 
#' @examples
#' # Mocking up a gene expression + ADT dataset:
#' library(scater)
#' exprs_sce <- mockSCE()
#' exprs_sce <- logNormCounts(exprs_sce)
#' exprs_sce <- runPCA(exprs_sce)
#' 
#' adt_sce <- mockSCE(ngenes=20) 
#' adt_sce <- logNormCounts(adt_sce)
#' altExp(exprs_sce, "ADT") <- adt_sce
#'
#' # Pretend we have three batches, and perform batch correction:
#' batch <- sample(1:3, ncol(exprs_sce), replace=TRUE)
#' corrected <- multiModalMNN(exprs_sce, batch=batch, which="ADT")
#'
#' # Pass arguments to, e.g., use a subset of features for the RNA,
#' # turn off the PCA for the ADTs:
#' corrected2 <- multiModalMNN(exprs_sce, batch=batch, which="ADT",
#'     pca.main.args=list(subset.row=1:500), 
#'     pca.alt.args=list(ADT=list(d=NA)))
#'
#' @export
#' @importFrom BiocNeighbors KmknnParam
#' @importFrom batchelor applyMultiSCE reducedMNN
#' @importFrom scuttle .unpackLists .bpNotSharedOrUp
#' @importFrom BiocParallel SerialParam bpstart bpstop
multiModalMNN <- function(..., batch=NULL, which=NULL, rescale.k=50, 
    pca.common.args=list(), pca.main.args=list(), pca.alt.args=list(),
    mnn.args=list(), BNPARAM=KmknnParam(), BPPARAM=SerialParam())
{
    if (.bpNotSharedOrUp(BPPARAM)) {
        bpstart(BPPARAM)
        on.exit(bpstop(BPPARAM), add=TRUE)
    }

    batches <- .unpackLists(...)
    for (x in seq_along(batches)) {
        if (!is(batches[[x]], "SingleCellExperiment")) {
            stop("entry ", x, " in '...' is not a SingleCellExperiment")
        }
    }

    all.pca <- applyMultiSCE(batches, FUN=.run_pca, WHICH=which, SIMPLIFY=FALSE,
        COMMON.ARGS=c(pca.common.args, list(batch=batch)),
        MAIN.ARGS=pca.main.args, ALT.ARGS=pca.alt.args)

    # Combining to rescale:
    combined <- all.pca
    if (length(batches) > 1) {
        for (i in seq_along(all.pca)) {
            combined[[i]] <- do.call(rbind, as.list(combined[[i]]))
        }
    } else {
        combined <- lapply(combined, "[[", i=1)
    }

    rescaled <- rescaleByNeighbors(combined, k=rescale.k, BNPARAM=BNPARAM, BPPARAM=BPPARAM)

    # Splitting it back for the MNN step.
    if (length(batches) > 1) {
        last <- 0L
        split.rescaled <- batches
        for (i in seq_along(batches)) {
            N <- ncol(batches[[i]])
            idx <- last + seq_len(N)
            split.rescaled[[i]] <- rescaled[idx,,drop=FALSE]
            last <- last + N
        }
        rescaled <- split.rescaled
    } else {
        rescaled <- list(rescaled)
    }

    corrected <- do.call(reducedMNN, c(rescaled, mnn.args, list(batch=batch, BNPARAM=BNPARAM, BPPARAM=BPPARAM)))
    .organize_outputs(all.pca, corrected)
}

#' @importFrom scuttle normalizeCounts
#' @importFrom batchelor cosineNorm multiBatchPCA
#' @importFrom SummarizedExperiment assay
#' @importFrom BiocParallel SerialParam
#' @importFrom BiocSingular IrlbaParam
.run_pca <- function(..., batch, cos.norm=TRUE, 
    d=50, 
    subset.row=NULL, 
    weights=NULL, 
    correct.all=FALSE, 
    get.variance=FALSE,
    assay.type="logcounts", 
    BSPARAM=IrlbaParam(),
    deferred=TRUE,
    BPPARAM=SerialParam()) 
{
    batches <- list(...)
    batches <- lapply(batches, assay, i=assay.type)

    if (cos.norm) {
        all.l2s <- lapply(batches, FUN=cosineNorm, mode="l2norm", subset.row=subset.row, BPPARAM=BPPARAM)
        for (i in seq_along(batches)) {
            l2 <- pmax(1e-8, all.l2s[[i]]) # protect against zero-L2.
            batches[[i]] <- normalizeCounts(batches[[i]], size.factors=l2, center.size.factors=FALSE, log=FALSE)
        }
    }

    do.call(multiBatchPCA, c(batches, 
        list(
            batch=batch, 
            d=d,
            subset.row=subset.row, 
            weights=weights,
            get.all.genes=correct.all,
            get.variance=get.variance,
            preserve.single=TRUE,
            BSPARAM=BSPARAM,
            deferred=deferred,
            BPPARAM=BPPARAM
        )
    ))
}

#' @importFrom S4Vectors metadata
#' @importFrom batchelor convertPCsToSCE
#' @importFrom SingleCellExperiment simplifyToSCE reducedDim<-
.organize_outputs <- function(all.pcs, corrected) {
    final <- vector("list", length(all.pcs))
    names(final) <- names(all.pcs)
    counter <- 0L

    for (i in seq_along(all.pcs)) {
        jump <- ncol(all.pcs[[i]][[1]])
        keep <- counter + seq_len(jump)
        tmp <- corrected
        tmp$corrected <- tmp$corrected[,keep,drop=FALSE]

        pc.mat <- all.pcs[[i]]
        pc.extras <- metadata(pc.mat)
        final[[i]] <- convertPCsToSCE(tmp, pc.extras, dimred.name=NULL)

        counter <- counter + jump
    }

    output <- simplifyToSCE(final, warn.level=3)
    reducedDim(output, "corrected") <- corrected$corrected
    output
}
