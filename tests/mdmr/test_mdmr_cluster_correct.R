library(connectir)
library(testthat)
library(stringr)

context("MDMR - Cluster Correct")

create_data <- function(nrow, ncol) {
    xm <- matrix(runif(nrow*ncol, min=0.01, max=10), nrow, ncol)
    xbm <- as.big.matrix(xm, type="double", shared=TRUE)
    return(xbm)
}

test_that("can get one column in big matrix", {
    bm <- create_data(10, 4)
    ref <- bm[,2]
    comp <- bm.get_col(bm, 2)
    expect_that(ref, is_equivalent_to(comp))
})

test_that("assign one column in big matrix with data", {
    bm <- create_data(10, 4)
    m <- as.matrix(bm)
    
    col <- rnorm(10)
    m[,2] <- col
    bm.save_col(col, bm, 2)
    
    ref <- m[,2]
    comp <- bm.get_col(bm, 2)
    expect_that(ref, is_equivalent_to(comp))
})

test_that("get p-values for each permutation", {
    fstats <- create_data(1000, 10)
    
    # Reference
    nperms <- nrow(fstats)
    nvoxs <- ncol(fstats)
    ref <- sapply(1:nvoxs, function(i) {
        ps <- (nperms - rank(fstats[,i], ties.method="max"))/nperms
        ps
    })
    
    # Comparison
    comp <- clust_mdmr.get_pvals_for_each_perm(fstats)
    
    expect_that(ref, is_equivalent_to(comp[,]))
})

test_that("get cluster details (mass, size, mass/size)", {
    nperms <- 1000
    nvoxs <- 10
    voxdim <- c(3,3,2)
    vox.thresh <- 0.5
    mask <- c(F, T, F, T, T, T, F, T, F, F, T, F, T, T, T, F, T, F)
    
    # Get p-values for each permutation
    fstats <- create_data(nperms, nvoxs)
    pvals <- clust_mdmr.get_pvals_for_each_perm(fstats)
    
    # Constrain permutations to cluster
    nperms <- 10
    pvals <- deepcopy(pvals, rows=1:nperms)
    

    # Reference
    ref <- laply(1:nperms, function(i) {
        ct <- cluster.table(1-pvals[i,], 1-vox.thresh, voxdim, mask)
        c(size=ct$max.size, mass=ct$max.mass, rel=ct$max.rel)
    })
    
    # Comparison
    comp <- clust_mdmr.values(pvals, mask, voxdim, vox.thresh=vox.thresh)
    
    expect_that(ref, is_equivalent_to(comp))
})


test_that("get cluster p-values for size", {
    nperms <- 1000
    nvoxs <- 10
    voxdim <- c(3,3,2)
    vox.thresh <- 0.5
    mask <- c(F, T, F, T, T, T, F, T, F, F, T, F, T, T, T, F, T, F)
    
    # Get p-values for each permutation
    fstats <- create_data(nperms, nvoxs)
    pvals <- clust_mdmr.get_pvals_for_each_perm(fstats)
    
    # Constrain permutations to cluster
    nperms <- 10
    pvals <- deepcopy(pvals, rows=1:nperms)
    
    # Get reference details & cluster values
    ref_details <- cluster.table(1-pvals[1,], 1-vox.thresh, voxdim, mask)
    clust_vals <- clust_mdmr.values(pvals, mask, voxdim, vox.thresh=vox.thresh)
    
    # Reference
    ref <- sapply(ref_details$size, function(x) sum(clust_vals[,1]>=x)/nperms)
    
    # Comparison
    comp <- clust_mdmr.pvalues(clust_vals, ref_details, "size")
    
    expect_that(ref, is_equivalent_to(comp$clust.pvals))
})


test_that("get cluster p-values for mass", {
    nperms <- 1000
    nvoxs <- 10
    voxdim <- c(3,3,2)
    vox.thresh <- 0.5
    mask <- c(F, T, F, T, T, T, F, T, F, F, T, F, T, T, T, F, T, F)
    
    # Get p-values for each permutation
    fstats <- create_data(nperms, nvoxs)
    pvals <- clust_mdmr.get_pvals_for_each_perm(fstats)
    
    # Constrain permutations to cluster
    nperms <- 10
    pvals <- deepcopy(pvals, rows=1:nperms)
    
    # Get reference details & cluster values
    ref_details <- cluster.table(1-pvals[1,], 1-vox.thresh, voxdim, mask)
    clust_vals <- clust_mdmr.values(pvals, mask, voxdim, vox.thresh=vox.thresh)
    
    # Reference
    ref <- sapply(ref_details$mass, function(x) sum(clust_vals[,2]>=x)/nperms)
    
    # Comparison
    comp <- clust_mdmr.pvalues(clust_vals, ref_details, "mass")
    
    expect_that(ref, is_equivalent_to(comp$clust.pvals))
})


test_that("find the significant clusters", {
    nperms <- 1000
    nvoxs <- 10
    voxdim <- c(3,3,2)
    vox.thresh <- 0.5
    clust.thresh <- 0.6
    mask <- c(F, T, F, T, T, T, F, T, F, F, T, F, T, T, T, F, T, F)
    
    # Get p-values for each permutation
    fstats <- create_data(nperms, nvoxs)
    pvals <- clust_mdmr.get_pvals_for_each_perm(fstats)
    
    # Constrain permutations to cluster
    nperms <- 10
    pvals <- deepcopy(pvals, rows=1:nperms)
    
    # Get reference details & cluster values
    clust_vals <- clust_mdmr.values(pvals, mask, voxdim, vox.thresh=vox.thresh)
    ref_details <- cluster.table(1-pvals[1,], 1-vox.thresh, voxdim, mask)
    ref_details <- clust_mdmr.pvalues(clust_vals, ref_details, "mass")
    
    # Reference
    clust <- ref_details$clust[mask]
    w.clusts <- which(rev(ref_details$clust.pvals<clust.thresh))
    clust.nums <- vector("numeric", length(clust))
    for (i in 1:length(w.clusts))
        clust.nums[clust==w.clusts[i]] <- i
    ref <- clust.nums
    
    # Comparison
    comp <- clust_mdmr.clusterize(ref_details, mask, clust.thresh)
    
    expect_that(ref, is_equivalent_to(comp))
})

test_that("cluster correct 'MDMR' results", {
    # Setup
    nperms <- 1000
    nvoxs <- 10
    hdr <- list(dim=c(3,3,2))
    vox.thresh <- 0.5
    clust.thresh <- 0.6
    mask <- c(F, T, F, T, T, T, F, T, F, F, T, F, T, T, T, F, T, F)
    
    # Pseudo-F Statistics
    fstats <- create_data(nperms, nvoxs)
    
    # Reference Clusters
    pvals <- clust_mdmr.get_pvals_for_each_perm(fstats)
    clust_vals <- clust_mdmr.values(pvals, mask, voxdim, vox.thresh=vox.thresh)
    ref_details <- cluster.table(1-pvals[1,], 1-vox.thresh, voxdim, mask)
    ref_details <- clust_mdmr.pvalues(clust_vals, ref_details, "mass")    
    ref <- clust_mdmr.clusterize(ref_details, mask, clust.thresh)
    
    # Comparison Clusters
    comp <- clust_mdmr.correct(mask, hdr, fstats, vox.thresh=vox.thresh, 
                        clust.thresh=clust.thresh, clust.type="mass")
    
    expect_that(ref, is_equivalent_to(comp))
})

test_that("qlm residuals works?", {
    dat <- create_data(k=2)
    
    fit <- lm(dat$y ~ dat$X - 1)
    residuals <- qlm_residuals(as.big.matrix(dat$y), as.big.matrix(dat$X))
    
    ref <- as.matrix(fit$residuals)
    comp <- as.matrix(residuals)
    expect_that(ref, is_equivalent_to(comp))
})
