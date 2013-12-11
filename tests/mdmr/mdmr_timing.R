# On rocky average was 38 secs and this was with accidental usage of multi-threading
# On gelert average was 30 secs and this was with single threading


###
# ROCKY
###

library(connectir)
library(testthat)
library(stringr)

create_init_data <- function(n) {
    formula <- ~ group*continuous
    model <- data.frame(group = rep(c("A","B","C","D"), each=n/4), 
                        continuous = rnorm(n))
    factor.names <- c("group", "continuous", "group:continuous")
    
    return(list(
        formula = formula, 
        model = model, 
        factor.names = factor.names
    ))
}

create_rhs <- function(n=100) {
    dat <- create_init_data(n)
    rhs <- mdmr_model.rhs(dat$formula, dat$model)
    return(rhs)
}

create_Gs <- function(n, nvoxs, outdir) {
    I <- diag(n)
    ones <- matrix(1, nrow=n)
    C <- I - ones %*% t(ones)/n
    Gs.mat <- sapply(1:nvoxs, function(vi) {
        smat <- matrix(rnorm(nvoxs*n), nvoxs, n)
        dmat <- 1 - cor(smat)
        A <- -(dmat^2)/2
        G <- C %*% A %*% C
        as.vector(G)
    })
    as.big.matrix(Gs.mat, backingpath=outdir, backingfile="gower.bin", descriptorfile="gower.desc")
}

nobs <- 100; nvoxs <- 12; nperms <- 100
dat <- create_init_data(nobs)
formula <- dat$formula
model <- dat$model
outdir <- tempdir()
Gs <- create_Gs(nobs, nvoxs, outdir)

reg.timings <- lapply(1:10, function(i) {
    duration <- system.time(
                    mdmr(Gs, formula, model, 
                         nperms=nperms, save.fperms=TRUE, G.path=outdir, 
                         blocksize=10, superblocksize=2, fperms.path=outdir)
                )
    rm_files <- list.files(pattern="fperms")
    file.remove(rm_files)
    duration
})




###
# GELERT
###

library(connectir)
library(testthat)
library(stringr)
library(Rsge)

create_init_data <- function(n) {
    formula <- ~ group*continuous
    model <- data.frame(group = rep(c("A","B","C","D"), each=n/4), 
                        continuous = rnorm(n))
    factor.names <- c("group", "continuous", "group:continuous")
    
    return(list(
        formula = formula, 
        model = model, 
        factor.names = factor.names
    ))
}

create_rhs <- function(n=100) {
    dat <- create_init_data(n)
    rhs <- mdmr_model.rhs(dat$formula, dat$model)
    return(rhs)
}

create_Gs <- function(n, nvoxs, outdir) {
    I <- diag(n)
    ones <- matrix(1, nrow=n)
    C <- I - ones %*% t(ones)/n
    Gs.mat <- sapply(1:nvoxs, function(vi) {
        smat <- matrix(rnorm(nvoxs*n), nvoxs, n)
        dmat <- 1 - cor(smat)
        A <- -(dmat^2)/2
        G <- C %*% A %*% C
        as.vector(G)
    })
    as.big.matrix(Gs.mat, backingpath=outdir, backingfile="gower.bin", descriptorfile="gower.desc")
}

nobs <- 100; nvoxs <- 12; nperms <- 100
dat <- create_init_data(nobs)
formula <- dat$formula
model <- dat$model
outdir <- file.path(getwd(), "tmp")
dir.create(outdir)
Gs <- create_Gs(nobs, nvoxs, outdir)

sge.timings <- lapply(1:10, function(i) {
    sge.info <- list(njobs=4, nforks=1, nthreads=1, ignore.proc.error=T)
    duration <- system.time(
                    mdmr(Gs, formula, model, 
                         nperms=nperms, save.fperms=TRUE, 
                         G.path=outdir, fperms.path=outdir, 
                         blocksize=10, superblocksize=2, sge.info=sge.info)
                 )
     rm_files <- list.files(outdir, pattern="fperms", full.names=T)
     file.remove(rm_files)
    duration
})
