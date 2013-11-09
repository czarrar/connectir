library(connectir)
library(testthat)
library(stringr)

context("correlations with big_cor")

# TODO: make non x %*% x tests!

create_data <- function(nr=100, nc=5, transpose=FALSE, scale=TRUE) {
    xm <- scale(matrix(rnorm(nr*nc), nr, nc), scale=scale)
    if (transpose)
        xm <- t(xm)
    xbm <- as.big.matrix(xm)
    return(list(m=xm, bm=xbm))
}

test_that("self correlation works", {
    dat <- create_data()
    ref <- cor(dat$m)
    comp <- big_cor(dat$bm)
    expect_that(ref, is_equivalent_to(as.matrix(comp)))
})

test_that("self inverse covariance works", {
    dat <- create_data(scale=FALSE)
    ref <- icov(dat$m)
    comp <- big_cor(dat$bm, glasso=T)
    expect_that(ref, is_equivalent_to(as.matrix(comp)))
})

test_that("one can pre-specify output matrix", {
    dat <- create_data()
    ref <- cor(dat$m)
    cmat <- big.matrix(ncol(dat$bm), ncol(dat$bm), shared=FALSE)
    big_cor(dat$bm, z=cmat)
    expect_that(ref, is_equivalent_to(as.matrix(cmat)))
})

test_that("self correlation with part of input works", {
    dat <- create_data(nc=10)
    si <- 5; ei <- 10
    ref <- cor(dat$m[,si:ei])
    comp <- big_cor(dat$bm, x_firstCol=5, x_lastCol=10)
    expect_that(ref, is_equivalent_to(as.matrix(comp)))
})

test_that("self correlation with part of output works", {
    dat <- create_data(nc=10)
    ref <- matrix(0, 10, 20)
    ref[,10:19] <- cor(dat$m)
    cmat <- big.matrix(10, 20, init=0, type="double", shared=FALSE)
    comp <- big_cor(dat$bm, z=cmat, z_firstCol=10, z_lastCol=19)
    expect_that(ref, is_equivalent_to(as.matrix(comp)))
})

test_that("error raised if output has incorrect number of rows", {
    dat <- create_data(nc=10)
    cmat <- big.matrix(20, 10, type="double")
    expect_that(big_cor(dat$bm, z=cmat), throws_error("incorrect number of rows"))
})

test_that("error raised if output has incorrect number of columns", {
    dat <- create_data(nc=10)
    cmat <- big.matrix(10, 20, type="double")
    expect_that(big_cor(dat$bm, z=cmat), throws_error("incorrect number of columns"))
})

test_that("self correlation by rows works", {
    dat <- create_data(nr=100, nc=5, transpose=TRUE)
    ref <- cor(t(dat$m))
    comp <- big_cor(dat$bm, byrows=TRUE)
    expect_that(ref, is_equivalent_to(as.matrix(comp)))
})

test_that("self correlation by rows with part of output works", {
    dat <- create_data(nr=100, nc=5, transpose=TRUE)
    ref <- matrix(0, 5, 20)
    ref[,5:9] <- cor(t(dat$m))
    cmat <- big.matrix(5, 20, init=0, type="double", shared=FALSE)
    comp <- big_cor(dat$bm, z=cmat, byrows=TRUE, z_firstCol=5, z_lastCol=9)
    expect_that(ref, is_equivalent_to(as.matrix(comp)))
})

test_that("wrapper vbca2 works", {
    dat <- create_data(nc=10)
    ref <- matrix(0, 10, 20)
    ref[,10:19] <- cor(dat$m)
    cmat <- big.matrix(10, 20, init=0, type="double", shared=FALSE)
    comp <- vbca2(dat$bm, outmat=cmat, outcols=c(10,19))
    expect_that(ref, is_equivalent_to(as.matrix(comp)))
})

test_that("wrapper vbca2 works with ztransform", {
    dat <- create_data(nc=5)
    ref <- atanh(cor(dat$m))
    comp <- vbca2(dat$bm, ztransform=TRUE)
    diag(ref) <- 0
    diag(comp) <- 0
    expect_that(ref, is_equivalent_to(as.matrix(comp)))
})

test_that("wrapper vbca2 works with ztransform when taking a subsample", {
    dat <- create_data(nc=5)
    ref <- atanh(cor(dat$m[,1:3], dat$m))
    comp <- vbca2(dat$bm, c(1,3), ztransform=TRUE)
    diag(ref) <- 0
    diag(comp) <- 0
    expect_that(ref, is_equivalent_to(as.matrix(comp)))
})

test_that("wrapper vbca3 works", {
    dat <- create_data(nc=10)
    ref <- matrix(0, 10, 20)
    ref[,10:19] <- cor(dat$m)
    cmat <- big.matrix(10, 20, init=0, type="double", shared=FALSE)
    comp <- vbca3(dat$bm, outmat=cmat, outcols=c(10,19))
    expect_that(ref, is_equivalent_to(as.matrix(comp)))
})

test_that("wrapper vbca3 works with 2 different inputs", {
    dat1 <- create_data(nc=10)
    dat2 <- create_data(nc=20)
    ref <- cor(dat1$m, dat2$m)
    cmat <- big.matrix(10, 20, init=0, type="double", shared=FALSE)
    comp <- vbca3(inmat1=dat1$bm, inmat2=dat2$bm)
    expect_that(ref, is_equivalent_to(as.matrix(comp)))
})

test_that("wrapper vbca3 works with ztransform", {
    dat <- create_data(nc=5)
    ref <- atanh(cor(dat$m))
    comp <- vbca3(dat$bm, ztransform=TRUE)
    diag(ref) <- 0
    diag(comp) <- 0
    expect_that(ref, is_equivalent_to(as.matrix(comp)))
})

test_that("wrapper vbca3 works with ztransform when taking a subsample", {
    dat <- create_data(nc=5)
    ref <- atanh(cor(dat$m[,1:3], dat$m))
    comp <- vbca3(dat$bm, c(1,3), ztransform=TRUE)
    diag(ref) <- 0
    diag(comp) <- 0
    expect_that(ref, is_equivalent_to(as.matrix(comp)))
})

