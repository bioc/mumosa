# Tests the unsplitAltExps function.
# library(testthat); library(mumosa); source("test-unsplit-altexp.R")

library(scuttle)
sce <- mockSCE()
se0 <- mockSCE(ngenes=20)
se <- as(se0, "SummarizedExperiment")
dimnames(se) <- dimnames(se0)
altExp(sce, "stuff") <- se

is.first <- seq_len(nrow(sce))
is.spike <- nrow(sce) + seq_len(nrow(altExp(sce)))
is.stuff <- nrow(sce) + nrow(altExp(sce)) + + seq_len(nrow(altExp(sce, 2)))

expected_names <- c(
    rownames(sce), 
    paste0("Spikes.", rownames(altExp(sce))),
    paste0("stuff.", rownames(altExp(sce, 2)))
)

test_that("unsplitAltExps works as expected for assays", {
    output <- unsplitAltExps(sce, prefix.row=FALSE)
    expect_identical(as.matrix(assay(output)), rbind(assay(sce), assay(altExp(sce)), assay(se)))

    output <- unsplitAltExps(sce)
    expect_identical(expected_names, rownames(output))

    # A trickier situation with an assay only available in one of the matrices.
    alt <- logNormCounts(sce)
    output <- unsplitAltExps(alt)
    expect_identical(as.matrix(logcounts(alt)), as.matrix(logcounts(output)[is.first,]))
    expect_true(all(is.na(logcounts(output)[is.spike,])))
    expect_true(all(is.na(logcounts(output)[is.stuff,])))

    alt <- sce
    altExp(alt) <- logNormCounts(altExp(alt))
    output <- unsplitAltExps(alt, prefix.rows=FALSE)
    expect_true(all(is.na(logcounts(output)[is.first,])))
    expect_identical(as.matrix(logcounts(altExp(alt))), as.matrix(logcounts(output)[is.spike,]))
    expect_true(all(is.na(logcounts(output)[is.stuff,])))

    # Counterpart is present in all assays.
    sce <- logNormCounts(sce, use.altexps=TRUE)
    output <- unsplitAltExps(sce, prefix.row=FALSE)
    expect_identical(as.matrix(logcounts(output)), rbind(logcounts(sce), logcounts(altExp(sce)), assay(altExp(sce, 2), "logcounts")))
})

test_that("rowRanges merges work as expected", {
    replacement <- rep(List(GRanges("chrA:1-100")), nrow(sce))
    names(replacement) <- rownames(sce)
    rowRanges(sce) <- replacement
    output <- unsplitAltExps(sce)

    len <- lengths(rowRanges(output))
    expect_identical(expected_names, names(len))
    expect_identical(unname(len), c(rep(1L, nrow(sce)), integer(nrow(altExp(sce))), integer(nrow(se))))

    # What happens if one of them is a GRangesList, and another is a GRanges?
    replacement <- GRanges(rep("chrB:1-100", nrow(altExp(sce))))
    names(replacement) <- rownames(altExp(sce))
    rowRanges(altExp(sce)) <- replacement
    suppressWarnings(output <- unsplitAltExps(sce))

    expect_s4_class(rowRanges(output), "GRangesList")
    is.A <- any(seqnames(rowRanges(output))=="chrA")
    expect_true(all(is.A[is.first]))
    expect_true(!any(is.A[-is.first]))
    is.B <- any(seqnames(rowRanges(output))=="chrB")
    expect_true(all(is.B[is.spike]))
    expect_true(!any(is.B[-is.spike]))
    expect_true(all(lengths(rowRanges(output))[is.stuff]==0))

    # What happens if all of them are GRanges?
    replacement <- GRanges(rep("chrA:1-100", nrow(sce)))
    names(replacement) <- rownames(sce)
    rowRanges(sce) <- replacement

    suppressWarnings(output <- unsplitAltExps(sce))
    expect_s4_class(rowRanges(output), "GRanges")
    expect_identical(as.character(seqnames(output)), 
        rep(c("chrA", "chrB", "unknown"), c(nrow(sce), nrow(altExp(sce)), nrow(altExp(sce, 2)))))

    # How are metadata fields handled?
    rowData(sce)$thingy <- 1
    rowData(altExp(sce))$blah <- "A"
    rowData(altExp(sce, 2))$whee <- TRUE

    suppressWarnings(output <- unsplitAltExps(sce))
    expect_type(rowData(output)$thingy, "double")
    expect_type(rowData(output)$blah, "character")
    expect_type(rowData(output)$whee, "logical")
})

test_that("colData merges work as expected", {
    output <- unsplitAltExps(sce)
    expect_identical(output$Cell_Cycle, sce$Cell_Cycle)
    expect_identical(output$stuff.Cell_Cycle, altExp(sce, 2)$Cell_Cycle)

    output <- unsplitAltExps(sce, prefix.cols=FALSE)
    expect_identical(sum(colnames(colData(output))=="Cell_Cycle"), 2L)
})

test_that("reducedDim merges work as expected", {
    reducedDim(sce, "PCA") <- matrix(rnorm(ncol(sce)*2), ncol=2)
    output <- unsplitAltExps(sce)
    expect_identical(reducedDims(sce), reducedDims(output))

    reducedDim(altExp(sce), "PCA") <- matrix(rnorm(ncol(sce)*2), ncol=2)
    output <- unsplitAltExps(sce)
    expect_identical(reducedDimNames(output), c("PCA", "Spikes.PCA"))
    expect_identical(reducedDim(output, "Spikes.PCA"), reducedDim(altExp(sce)))

    output <- unsplitAltExps(sce, prefix.cols=FALSE)
    expect_identical(reducedDimNames(output), c("PCA", "PCA"))
    expect_identical(reducedDim(output, 2), reducedDim(altExp(sce)))
})
