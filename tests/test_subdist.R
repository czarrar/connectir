library(connectir)
library(testthat)
library(stringr)
library(glasso)

context("Subject Distances")

create_data <- function(n=10, ncols=5) {
    xm <- matrix(runif((n^2)*ncols, min=0, max=2), n^2, ncols)
    xbm <- as.big.matrix(xm, type="double", shared=FALSE)
    return(list(m=xm, bm=xbm))
}

create_scale_data <- function(n=100, ncols=5) {
    xm <- matrix(runif((n^2)*ncols, min=0, max=2), n^2, ncols)
    xbm <- as.big.matrix(xm, type="double", shared=FALSE)
    return(list(m=xm, bm=xbm))
}

create_many_big_matrices <- function(nsubs=10, xdim=20, ydim=10, ...) {
    lapply(1:nsubs, function(x) as.big.matrix(scale(matrix(rnorm(xdim*ydim), xdim, ydim)), ...))
}

# G <- -0.5 * dmat^2 %*% (I - ones %*% t(ones)/n)
test_that("creation of gower matrix works", {
    dat <- create_data()
    
    n <- sqrt(nrow(dat$m))
    I <- diag(n)
    ones <- matrix(1, n, 1)
    adj <- I - tcrossprod(ones)/n
    ref <- apply(dat$m, 2, function(x) {
        dmat <- matrix(x, n, n)
        A <- -1*(dmat^2)/2
        as.vector(A %*% adj)
    })
    
    G <- big.matrix(nrow(dat$bm), ncol(dat$bm), type="double", shared=FALSE)
    nc <- ncol(dat$bm)
    comp <- .Call("big_gower", dat$bm, G, as.double(1), as.double(nc), 
                  as.double(1), as.double(nc), PACKAGE="connectir")
    expect_that(ref, is_equivalent_to(as.matrix(G)))
})

test_that("combining subject cormaps works", {
    firstSeed <- 51; lastSeed <- 60
    seedinds <- firstSeed:lastSeed; nseeds <- length(seedinds)
    voxinds <- 1:100; nvoxs <- length(voxinds)
    nsubs <- 20
    
    seedi <- 1
    
    # CorMaps
    funcs <- create_many_big_matrices(nsubs, 100, nvoxs, type="double", shared=FALSE)
    cormaps <- vbca_batch2(funcs, c(firstSeed, lastSeed), ztransform=TRUE, shared=FALSE)
    
    # Ref
    seedMapsA <- scale(sapply(cormaps, function(x) x[seedi,-voxinds[seedinds[seedi]]]))
    
    # Comp
    seedMapsB <- big.matrix(nvoxs-1, nsubs, type="double", shared=FALSE)
    tmp <- .Call("subdist_combine_and_scale_submaps", cormaps, 
          as.double(seedi), as.double(voxinds[-seedinds[seedi]]), seedMapsB, 
          PACKAGE="connectir")
    
    expect_that(seedMapsA, is_equivalent_to(as.matrix(seedMapsB)))
})

test_that("pearson distance works", {
    nvoxs <- 100; nsubs <- 10
    seedCorMaps <- as.big.matrix(scale(matrix(rnorm(nvoxs*nsubs), nvoxs, nsubs)))
    
    # Ref
    dmatA <- 1 - cor(seedCorMaps[,])
    
    # Comp
    dmatB <- big.matrix(nsubs^2, 1, type="double", shared=FALSE)
    .Call("subdist_pearson_distance", seedCorMaps, dmatB, as.double(1), FALSE, 
          PACKAGE="connectir")
    
    expect_that(as.vector(dmatA), is_equivalent_to(dmatB[,]))
})

test_that("sparse inverse covariance distance works", {
    nvoxs <- 100; nsubs <- 10
    seedCorMaps <- as.big.matrix(matrix(rnorm(nvoxs*nsubs), nvoxs, nsubs))
	
    # Ref
	tmp <- scale(seedCorMaps[,], center=TRUE, scale=FALSE)
    dmatA <- 1 - icov(tmp)
    
    # Comp
	center_fast(seedCorMaps, to.copy=FALSE, byrows=FALSE)
	oc <- big_cor(seedCorMaps, byrows=FALSE)
	dmatB <- 1 - norm_glasso(oc[,])
	
    expect_that(dmatA, equals(dmatB))
})

test_that("sparse inverse covariance distance works with transposed matrix", {
    nvoxs <- 100; nsubs <- 10
    seedCorMaps <- as.big.matrix(matrix(rnorm(nvoxs*nsubs), nsubs, nvoxs))
	
    # Ref
	tmp <- scale(t(seedCorMaps[,]), center=TRUE, scale=FALSE)
    dmatA <- 1 - icov(tmp)
    
    # Comp
	center_fast(seedCorMaps, to.copy=FALSE, byrows=TRUE)
	oc <- big_cor(seedCorMaps, byrows=TRUE)
	dmatB <- 1 - norm_glasso(t(oc[,]))
	
    expect_that(dmatA, equals(dmatB))
})

test_that("subdist worker works", {
    firstSeed <- 51; lastSeed <- 60
    seedinds <- firstSeed:lastSeed; nseeds <- length(seedinds)
    firstDist <- 1; lastDist <- 10
    distinds <- firstDist:lastDist; ndists <- length(distinds)
    voxinds <- 1:100; nvoxs <- length(voxinds)
    nsubs <- 20
    if (ndists != nseeds)
        stop("length mismatch")
    
    funcs <- create_many_big_matrices(nsubs, 100, nvoxs, type="double", shared=FALSE)
    
    # Ref
    cormaps <- vbca_batch2(funcs, c(firstSeed, lastSeed), ztransform=TRUE, shared=FALSE)
    dmatsA <- big.matrix(nsubs^2, ndists, type="double", shared=FALSE)
    for (i in 1:nseeds) {
        seedMaps <- scale(sapply(cormaps, function(x) x[i,-voxinds[seedinds[i]]]))
        dmatsA[,distinds[i]] <- as.vector(1 - cor(seedMaps))
    }
    
    # Comp
    dmatsB <- big.matrix(nsubs^2, ndists, type="double", shared=FALSE)
    compute_subdist_worker2(funcs, firstSeed, lastSeed, 
                            dmatsB, firstDist, lastDist, 
                            ztransform=TRUE, method="pearson", 
                            type="double", shared=FALSE)
    
    expect_that(dmatsA[,], equals(dmatsB[,]))
})

test_that("subdist worker with regression first works", {
    firstSeed <- 51; lastSeed <- 60
    seedinds <- firstSeed:lastSeed; nseeds <- length(seedinds)
    firstDist <- 1; lastDist <- 10
    distinds <- firstDist:lastDist; ndists <- length(distinds)
    voxinds <- 1:100; nvoxs <- length(voxinds)
    nsubs <- 20
    if (ndists != nseeds)
        stop("length mismatch")
    
    funcs <- create_many_big_matrices(nsubs, 100, nvoxs, type="double", shared=FALSE)
    
    design_mat <- cbind(rep(1,nsubs), sample(0:1, nsubs, replace=T), rnorm(nsubs))
    
    # Ref
    cormaps <- vbca_batch2(funcs, c(firstSeed, lastSeed), ztransform=TRUE, shared=FALSE)
    dmatsA <- big.matrix(nsubs^2, ndists, type="double", shared=FALSE)
    for (i in 1:nseeds) {
        seedMaps <- sapply(cormaps, function(x) x[i,-voxinds[seedinds[i]]])
        fit <- lm(t(seedMaps) ~ design_mat)
        seedMaps <- scale(t(fit$residuals))
        dmatsA[,distinds[i]] <- as.vector(1 - cor(seedMaps))
    }
    
    # Comp
    dmatsB <- big.matrix(nsubs^2, ndists, type="double", shared=FALSE)
    compute_subdist_worker2_regress(funcs, firstSeed, lastSeed, 
                                    dmatsB, firstDist, lastDist, 
                                    as.big.matrix(design_mat), 
                                    ztransform=TRUE, method="pearson", 
                                    type="double", shared=FALSE)
    
    expect_that(dmatsA[,], equals(dmatsB[,]))
})

test_that("subdist works", {
    firstSeed <- 51; lastSeed <- 70
    seedinds <- firstSeed:lastSeed; nseeds <- length(seedinds)
    firstDist <- 1; lastDist <- nseeds
    distinds <- firstDist:lastDist; ndists <- length(distinds)
    voxinds <- 1:100; nvoxs <- length(voxinds)
    nsubs <- 20
    
    funcs <- create_many_big_matrices(nsubs, 100, nvoxs, type="double", shared=FALSE)
    
    # Ref
    cormaps <- vbca_batch2(funcs, c(firstSeed, lastSeed), ztransform=TRUE, shared=FALSE)
    dmatsA <- big.matrix(nsubs^2, ndists, type="double", shared=FALSE)
    for (i in 1:nseeds) {
        seedMaps <- scale(sapply(cormaps, function(x) x[i,-voxinds[seedinds[i]]]))
        dmatsA[,distinds[i]] <- as.vector(1 - cor(seedMaps))
    }
    
    # Comp
    dmatsB <- big.matrix(nsubs^2, ndists, type="double", shared=FALSE)
    compute_subdist2(funcs, firstSeed, lastSeed, 
                    dmatsB, firstDist, lastDist, 
                    blocksize=5, ztransform=TRUE, method="pearson", 
                    type="double")
    
    expect_that(dmatsA[,], equals(dmatsB[,]))
})
