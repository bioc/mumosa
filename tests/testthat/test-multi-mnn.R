# This tests the multiModalMNN function.
# library(testthat); library(mumosa); source("test-multi-mnn.R")

set.seed(1000)
library(scater)
exprs_sce <- mockSCE()
exprs_sce <- logNormCounts(exprs_sce)
exprs_sce <- runPCA(exprs_sce)

adt_sce <- mockSCE(ngenes=20) 
adt_sce <- logNormCounts(adt_sce)
altExp(exprs_sce, "ADT") <- adt_sce

BSPARAM <- BiocSingular::ExactParam()

test_that("multiModalMNN recapitulates fastMNN", {
    batch <- sample(1:3, ncol(exprs_sce), replace=TRUE)

    corrected <- multiModalMNN(exprs_sce, batch=batch, which=character(0), common.args=list(BSPARAM=BSPARAM))
    ref <- batchelor::fastMNN(exprs_sce, batch=batch, BSPARAM=BSPARAM)

    ratio <- reducedDim(corrected)/reducedDim(ref)
    expect_true(mad(ratio) < 1e-8)
    expect_equal(as.matrix(assay(corrected)), as.matrix(assay(ref) * median(ratio)))

    # Trying with multiple experiments.
    copy <- exprs_sce
    altExp(copy) <- removeAltExps(exprs_sce)
    corrected2 <- multiModalMNN(copy, batch=batch, which="Spikes", common.args=list(BSPARAM=BSPARAM))

    expect_equal(reducedDim(corrected), reducedDim(corrected2)[,1:50])
    expect_equal(reducedDim(corrected), reducedDim(corrected2)[,50+1:50])
    expect_equal(assay(corrected), assay(corrected2))
    expect_equal(assay(corrected), assay(altExp(corrected2)))

    # Trying with experiments split across multiple objects.
    corrected3 <- multiModalMNN(
        exprs_sce[,batch==1], exprs_sce[,batch==2], exprs_sce[,batch==3],
        which=character(0), common.args=list(BSPARAM=BSPARAM)
    )

    ref <- reducedDim(corrected)[order(batch),]
    obs <- reducedDim(corrected3)
    ratio2 <- abs(ref/obs)
    expect_equal(mean(ratio2), 1)
    expect_true(sd(ratio2) < 1e-8)

    # And again, with multiple objects in multiple experiments.
    corrected4 <- multiModalMNN(
        copy[,batch==1], copy[,batch==2], copy[,batch==3],
        which="Spikes", common.args=list(BSPARAM=BSPARAM)
    )

    ref <- reducedDim(corrected2)[order(batch),]
    obs <- reducedDim(corrected4)
    ratio2 <- abs(ref/obs)
    expect_equal(mean(ratio2), 1)
    expect_true(sd(ratio2) < 1e-8)
})

test_that("multiModalMNN rescales correctly", {
    batch <- sample(1:3, ncol(exprs_sce), replace=TRUE)

    suppressWarnings(ref <- multiModalMNN(exprs_sce, batch=batch, which="ADT",
        common.args=list(BSPARAM=BSPARAM)))

    copy <- exprs_sce
    logcounts(altExp(copy,2 )) <- 10*logcounts(altExp(copy, 2))
    suppressWarnings(alt <- multiModalMNN(copy, batch=batch, which="ADT",
        common.args=list(BSPARAM=BSPARAM)))

    expect_equal(alt, ref)
})

test_that("multiModalMNN respects Experiment-specific arguments", {
    batch <- sample(1:3, ncol(exprs_sce), replace=TRUE)

    corrected <- multiModalMNN(exprs_sce, batch=batch, which="ADT",
        common.args=list(BSPARAM=BSPARAM), main.args=list(subset.row=1:100), 
        alt.args=list(ADT=list(d=NA)))

    ref <- multiModalMNN(exprs_sce[1:100,], batch=batch, which="ADT",
        common.args=list(BSPARAM=BSPARAM), 
        alt.args=list(ADT=list(d=NA)))

    expect_identical(ref, corrected)
    expect_identical(as.matrix(assay(altExp(ref))), tail(t(reducedDim(ref)), 20)) # expression values directly used in reduced dims.

    # Smart enough to correct all.
    corrected2 <- multiModalMNN(exprs_sce, batch=batch, which="ADT",
        common.args=list(BSPARAM=BSPARAM), 
        main.args=list(subset.row=1:100, correct.all=TRUE), 
        alt.args=list(ADT=list(d=NA)))

    expect_identical(reducedDim(corrected), reducedDim(corrected2))
    expect_identical(assay(corrected), assay(corrected2)[1:100,])
    expect_identical(nrow(corrected2), nrow(exprs_sce))
})

test_that("multiModalMNN throws errors in the right places", {
    batch <- sample(1:3, ncol(exprs_sce), replace=TRUE)
    expect_error(corrected <- multiModalMNN(exprs_sce, batch=batch, which="ADT", BSPARAM=BSPARAM), 'not a SingleCellExperiment')
})
