library(connectir)
library(testthat)
library(stringr)
library(bigmemory)

context("kendall")

# input for different tests
create_many_big_matrices <- function(nsubs=10, xdim=20, ydim=10, ...) {
    lapply(1:nsubs, function(x) as.big.matrix(scale(matrix(rnorm(xdim*ydim), xdim, ydim)), ...))
}

# simpler function
simple_kendall <- function(bigmats) {
    nsubs <- length(bigmats)
    nvoxs <- ncol(bigmats[[1]])
    cormats <- vbca_batch(bigmats, 1:nvoxs)
    coeffs <- sapply(1:nvoxs, function(i) {
        xmat <- sapply(cormats, function(x) x[i,-i])
        kendall_ref(xmat)
    })
    return(coeffs)
}
simple_kendall3 <- function(bigmats, bigmats2=NULL) {
    nsubs <- length(bigmats)
    nvoxs <- ncol(bigmats[[1]])
    cormats <- vbca_batch3(bigmats, c(1,nvoxs), bigmats2)
    coeffs <- sapply(1:nvoxs, function(i) {
        xmat <- sapply(cormats, function(x) x[i,])
        kendall_ref(xmat)
    })
    return(coeffs)
}

simple_kendall3_regress <- function(bigmats, design_mat, bigmats2=NULL) {
    nsubs <- length(bigmats)
    nvoxs <- ncol(bigmats[[1]])
    cormats <- vbca_batch3(bigmats, c(1,nvoxs), bigmats2)
    coeffs <- sapply(1:nvoxs, function(i) {
        xmat <- sapply(cormats, function(x) x[i,])
        fit <- lm(t(xmat) ~ design_mat)
        r_xmat <- t(fit$residuals) + (rowMeans(xmat) %*% matrix(1, 1, nsubs))
        kendall_ref(r_xmat)
    })
    return(coeffs)
}

# test
test_that("kendall works", {
    xm <- matrix(rnorm(200*100), 200, 100)
    xbm <- as.big.matrix(xm)
    ref <- kendall_ref(xm)
    comp <- .Call("kendall_worker", xbm)
    ns <- nrow(xbm); nr <- ncol(xbm)
    comp <- (12 * comp * (ns - 1))/(nr^2 * (ns^3 - ns))
    expect_that(ref, equals(comp))
})

test_that("seedmap in kendall voxelwise works", {
    xs <- create_many_big_matrices(nsubs=10, xdim=200, ydim=100)
    
    seeds <- as.double(1:100); voxs <- seeds
    comp <- big.matrix(100-1, 10, shared=FALSE, type="double")
    cormats <- vbca_batch2(xs, c(1,100))
    tmp <- .Call("voxelwise_kendall", cormats, comp, seeds, voxs)
    
    ref <- sapply(cormats, function(x) x[100,-100])
    
    # comp <- kendall(xs, blocksize=10)
    expect_that(ref, is_equivalent_to(comp[,]))
})

test_that("kendall voxelwise works #1", {
    xs <- create_many_big_matrices(nsubs=10, xdim=200, ydim=100)

    seeds <- as.double(1:100); voxs <- seeds
    seedMaps <- big.matrix(100-1, 10, shared=FALSE, type="double")
    cormats <- vbca_batch2(xs, c(1,100))
    comp <- .Call("voxelwise_kendall", cormats, seedMaps, seeds, voxs)

    ref <- simple_kendall(xs)
    
    expect_that(ref, is_equivalent_to(comp[,]))
})

test_that("kendall voxelwise works #2", {
    xs <- create_many_big_matrices(nsubs=10, xdim=200, ydim=100)

    ref <- simple_kendall(xs)
    comp <- kendall(xs, blocksize=10)
    
    expect_that(ref, equals(comp))
})

test_that("kendall3 voxelwise works #1", {
    xs <- create_many_big_matrices(nsubs=10, xdim=200, ydim=100)

    ref <- simple_kendall3(xs)
    comp <- kendall3(xs, blocksize=10)
    
    expect_that(ref, equals(comp))
})

test_that("kendall3 voxelwise works #2", {
    xs <- create_many_big_matrices(nsubs=10, xdim=200, ydim=100)
    xs2 <- create_many_big_matrices(nsubs=10, xdim=200, ydim=200)
    
    ref <- simple_kendall3(xs, xs2)
    comp <- kendall3(xs, blocksize=10, xs2)
    
    expect_that(ref, equals(comp))
})

test_that("kendall3 voxelwise with regression works #1", {
    xs <- create_many_big_matrices(nsubs=10, xdim=200, ydim=100)
    design_mat <- matrix(rnorm(10*2), 10, 2)
    design_mat <- cbind(rep(1, 10), design_mat)
    
    ref0 <- simple_kendall3(xs)
    ref <- simple_kendall3_regress(xs, design_mat)
    comp <- kendall3(xs, blocksize=10, 
                    design_mat=as.big.matrix(design_mat, type="double"))
    
    expect_that(ref, equals(comp))
})
