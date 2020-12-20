# This tests the computeCorrelations function.
# library(testthat); library(mumosa); source("test-compute-corr.R")

library(scuttle)

set.seed(10001)
test_that("computeCorrelations works correctly", {
    sce1 <- mockSCE()
    sce1 <- logNormCounts(sce1)
    sce2 <- mockSCE(ngenes=10) 
    sce2 <- logNormCounts(sce2)

    output <- computeCorrelations(sce1, sce2)
    expect_identical(nrow(output), nrow(sce1) * nrow(sce2))

    # Randomly select 100 pairs.
    chosen <- sample(nrow(output), 100)
    vec.rho <- vec.p <- numeric(length(chosen))
    for (i in seq_along(chosen)) {
        j <- chosen[i]
        ct <- cor.test(
            logcounts(sce1)[output$feature1[j],], 
            logcounts(sce2)[output$feature2[j],], 
            method="spearman",
            exact=FALSE
        )
        vec.rho[i] <- ct$estimate
        vec.p[i] <- ct$p.value
    }
    expect_identical(vec.rho, output$rho[chosen])
    expect_identical(vec.p, output$p.value[chosen])

    # Handles subset.cols correctly
    sub <- computeCorrelations(sce1, sce2, subset.cols=10:100)
    ref <- computeCorrelations(sce1[,10:100], sce2[,10:100])
    expect_identical(sub, ref)
})

set.seed(10002)
test_that("computeCorrelations works correctly by itself", {
    sce <- mockSCE(ngenes=20)
    sce <- logNormCounts(sce)

    ref <- computeCorrelations(sce, sce)
    ref <- ref[ref$feature1!=ref$feature2,]
    ref$FDR <- p.adjust(ref$p.value, method="BH")

    output <- computeCorrelations(sce, y=NULL)
    expect_identical(ref, output)

    # Works correctly with subset.cols.
    sub <- computeCorrelations(sce, y=NULL, subset.cols=10:100)
    ref <- computeCorrelations(sce[,10:100], y=NULL)
    expect_identical(sub, ref)
})

set.seed(10003)
test_that("computeCorrelations works correctly with blocking", {
    sce1 <- mockSCE()
    sce1 <- logNormCounts(sce1)
    sce2 <- mockSCE(ngenes=10) 
    sce2 <- logNormCounts(sce2)

    ref <- computeCorrelations(sce1, sce2)
    output <- computeCorrelations(sce1, sce2, block=rep(1, ncol(sce1)))
    expect_identical(ref, output)

    # Equiweighting preserves the results when blocks are unbalanced.
    block <- rep(1:2, each=ncol(sce1)/2)
    ref <- computeCorrelations(sce1, sce2, block=block)

    expanded <- c(seq_len(ncol(sce1)), seq_len(ncol(sce1)/2))
    out <- computeCorrelations(sce1[,expanded], y=sce2[,expanded], block=block[expanded])

    ref <- ref[,1:3]
    ref <- ref[do.call(order, as.list(ref)),]
    out <- out[,1:3]
    out <- out[do.call(order, as.list(out)),]
    expect_equal(ref, out)

    # Same results without equiweighting, for balanced blocks.
    weight <- computeCorrelations(sce1, y=sce2, block=block)
    noweight <- computeCorrelations(sce1, y=sce2, block=block, equiweight=FALSE)
    expect_equal(weight, noweight)

    # Not true of unweighted and unbalanced blocks.    
    block0 <- rep(1:2, c(50, ncol(sce1)-50))
    weight <- computeCorrelations(sce1, y=sce2, block=block0)
    noweight <- computeCorrelations(sce1, y=sce2, block=block0, equiweight=FALSE)
    expect_false(isTRUE(all.equal(weight, noweight)))

    # subset.cols works as expected.
    sub <- computeCorrelations(sce1, sce2, block=block, subset.cols=10:100)
    ref <- computeCorrelations(sce1[,10:100], sce2[,10:100], block=block[10:100])
    expect_identical(sub, ref)
})

library(DelayedArray)
set.seed(10004)
test_that("computeCorrelations handles the looping correctly", {
    sce1 <- mockSCE(ngenes=20)
    counts(sce1) <- DelayedArray(counts(sce1))
    sce1 <- logNormCounts(sce1)

    sce2 <- mockSCE(ngenes=25)
    counts(sce2) <- DelayedArray(counts(sce2))
    sce2 <- logNormCounts(sce1)

    ref1 <- computeCorrelations(sce1, y=NULL)
    ref2 <- computeCorrelations(sce1, sce2)

    old <- getAutoBlockSize()
    setAutoBlockSize(ncol(sce1)*8*2)

    out1 <- computeCorrelations(sce1, y=NULL)
    out2 <- computeCorrelations(sce1, sce2)

    setAutoBlockSize(old)

    expect_identical(out1, ref1)
    expect_identical(out2, ref2)
})

set.seed(10005)
test_that("computeCorrelations handles the names correctly", {
    sce1 <- mockSCE()
    sce1 <- logNormCounts(sce1)
    sce2 <- mockSCE(ngenes=10) 
    sce2 <- logNormCounts(sce2)

    ref <- computeCorrelations(sce1, sce2)
    num.only <- computeCorrelations(sce1, sce2, use.names=FALSE)
    expect_identical(ref$feature1, rownames(sce1)[num.only$feature1])
    expect_identical(ref$feature2, rownames(sce2)[num.only$feature2])

    rowData(sce1)$BLAH <- tolower(rownames(sce1))
    rowData(sce2)$BLAH <- tolower(rownames(sce2))
    use.rd <- computeCorrelations(sce1, sce2, use.names="BLAH")
    expect_identical(use.rd$feature1, tolower(ref$feature1))
    expect_identical(use.rd$feature2, tolower(ref$feature2))

    use.rd <- computeCorrelations(sce1, sce2, use.names=c(NA, "BLAH"))
    expect_identical(use.rd$feature1, ref$feature1)
    expect_identical(use.rd$feature2, tolower(ref$feature2))

    use.rd <- computeCorrelations(sce1, logcounts(sce2), use.names="BLAH")
    expect_identical(use.rd$feature1, tolower(ref$feature1))
    expect_identical(use.rd$feature2, ref$feature2)
})
