library(connectir)
library(testthat)
library(stringr)

context("MDMR")

create_data <- function(nr, nc) {
    xm <- matrix(runif(nr*nc, min=0, max=2), nr, nc)
    xbm <- as.big.matrix(xm, type="double", shared=FALSE)
    return(xbm)
}

# G <- -0.5 * dmat^2 %*% (I - ones %*% t(ones)/n)
test_that("C++ mdmr_worker works", {
    nvoxs <- 100
    nsubs <- 20
    nperms <- 5000
    nfactors <- 2

    # Inputs
    Gmat <- create_data(nsubs^2, nvoxs)
    H2mats <- lapply(1:nfactors, function(i) create_data(nsubs^2, nperms))
    IHmat <- create_data(nsubs^2, nperms)
    dfRes <- as.double(nsubs - nfactors)
    dfExp <- as.double(rep(1, nfactors))
    
    # Outputs
    rFperms <- lapply(1:nfactors, function(i) create_data(nperms, nvoxs))
    cFperms <- lapply(rFperms, deepcopy, shared=FALSE, type="double")
    rPmat <- big.matrix(nvoxs, nfactors, init=1, type="double", shared=FALSE)
    cPmat <- deepcopy(rPmat, shared=FALSE, type="double")
    
    # Reference
    ErrorVariance <- t(IHmat[,]) %*% Gmat[,];
    for (i in 1:nfactors) {
        ExplainedVariance <- t(H2mats[[i]][,]) %*% Gmat[,]
        Fmat <- (ExplainedVariance/ErrorVariance) * (dfRes/dfExp[i])
        realFs <- Fmat[1,]
        rPmat[,i] <- sapply(1:nvoxs, function(j) sum(Fmat[,j] >= realFs[j]))
        rFperms[[i]][,] <- Fmat
    }
    
    # Comparison
    .Call("mdmr_worker", Gmat, cFperms, cPmat, 
          H2mats, IHmat, dfRes, dfExp, PACKAGE="connectir")
    
    expect_that(cFperms[[1]][,], equals(rFperms[[1]][,]))
    expect_that(cFperms[[2]][,], equals(rFperms[[2]][,]))
    expect_that(cPmat[,], equals(rPmat[,]))
})

