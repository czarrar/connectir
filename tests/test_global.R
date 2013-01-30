library(connectir)
library(testthat)
library(stringr)

context("global connectivity")

# input for different tests
create_cormat <- function(xdim=200, ydim=100, shared=FALSE, ...) {
    tmp <- matrix(rnorm(xdim*ydim), xdim, ydim)
    m <- cor(tmp)
    bm <- as.big.matrix(m, shared=shared, ...)
    list(m=m, bm=bm)
}

create_data <- function(xdim=200, ydim=100, shared=FALSE, ...) {
    tmp <- scale(matrix(rnorm(xdim*ydim), xdim, ydim))
    m <- cor(tmp)
    bm <- as.big.matrix(tmp, shared=shared, ...)
    list(m=m, bm=bm)
}

# simpler function
global_ref <- function(cormat, thresh=0, threshType=0) {
    diag(cormat) <- NA
    
    if (threshType == 0) {
        gcor <- colMeans(cormat, na.rm=TRUE)
    } else if (threshType == 1) {
        gcor <- apply(cormat, 2, function(x) mean(x[x>thresh], na.rm=TRUE))
    } else if (threshType == 2) {
        gcor <- apply(cormat, 2, function(x) sum((x>thresh)*1, na.rm=TRUE))
    } else if (threshType == 3) {
        gcor <- apply(cormat, 2, function(x) mean(x[x<thresh], na.rm=TRUE))
    } else if (threshType == 4) {
        gcor <- apply(cormat, 2, function(x) sum((x<thresh)*1, na.rm=TRUE))
    } else {
        stop("threshType not recognized")
    }
    
    return(gcor)
}

# test rowsums
test_that("gcor_worker works with thresh-type of 0", {
    l <- create_cormat(ydim=100)
    ref <- global_ref(l$m)
    comp <- .Call("gcor_worker", l$bm, as.double(0), as.integer(0), 
                  as.double(1:100), as.double(1:100), package="connectir")
    expect_that(ref, equals(comp[,1]))
})

test_that("gcor_worker works with thresh-type of 1", {
    l <- create_cormat(ydim=100)
    ref <- global_ref(l$m, 0, 1)
    comp <- .Call("gcor_worker", l$bm, as.double(0), as.integer(1), 
                  as.double(1:100), as.double(1:100), package="connectir")
    expect_that(ref, equals(comp[,1]))
})

test_that("gcor_worker works with thresh-type of 2", {
    l <- create_cormat(ydim=100)
    ref <- global_ref(l$m, 0, 2)
    comp <- .Call("gcor_worker", l$bm, as.double(0), as.integer(2), 
                  as.double(1:100), as.double(1:100), package="connectir")
    expect_that(ref, equals(comp[,1]))
})

test_that("gcor_worker works with thresh-type of 3", {
    l <- create_cormat(ydim=100)
    ref <- global_ref(l$m, 0, 3)
    comp <- .Call("gcor_worker", l$bm, as.double(0), as.integer(3), 
                  as.double(1:100), as.double(1:100), package="connectir")
    expect_that(ref, equals(comp[,1]))
})

test_that("gcor_worker works with thresh-type of 4", {
    l <- create_cormat(ydim=100)
    ref <- global_ref(l$m, 0, 4)
    comp <- .Call("gcor_worker", l$bm, as.double(0), as.integer(4), 
                  as.double(1:100), as.double(1:100), package="connectir")
    expect_that(ref, equals(comp[,1]))
})

# test main function
test_that("gcor works", {
    l <- create_data(ydim=100)
    ref <- global_ref(l$m)
    comp <- gcor(l$bm, 10, verbosity=2)
    
    expect_that(ref, equals(comp))
})

test_that("gcor with ztransform works", {
    l <- create_data(ydim=100)
    ref <- tanh(global_ref(atanh(l$m)))
    comp <- gcor(l$bm, 10, verbosity=2, ztransform=TRUE)
    
    expect_that(ref, equals(comp))
})
