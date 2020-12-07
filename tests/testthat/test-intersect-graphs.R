# This tests the graph intersection functionality.
# library(testthat); library(mumosa); source("test-intersect-graphs.R")

library(scran)

mat1 <- matrix(rnorm(20000), ncol=20)
chosen <- 1:250
mat1[chosen,1] <- mat1[chosen,1] + 10
g1 <- buildSNNGraph(mat1, d=NA, transposed=TRUE)

# Pretending we have some other data for the same cells, e.g., ADT.
mat2 <- matrix(rnorm(20000), ncol=20)
chosen <- c(1:125, 251:375)
mat2[chosen,2] <- mat2[chosen,2] + 10
g2 <- buildSNNGraph(mat2, d=NA, transposed=TRUE)

test_that("intersectGraphs works correctly", {
    G1 <- g1[]
    nz1 <- Matrix::which(G1 > 0)
    G2 <- g2[]    
    nz2 <- Matrix::which(G2 > 0)

    G1[setdiff(nz2, nz1)] <- min(G1[nz1]) * 1e-6
    G2[setdiff(nz1, nz2)] <- min(G2[nz2]) * 1e-6
    prod <- G1 * G2

    gcom <- intersectGraphs(g1, g2)
    ref <- gcom[]
    expect_identical(ref, prod)
})
