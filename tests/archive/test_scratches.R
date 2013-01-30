library(connectir)
library(testthat)
library(stringr)

context("testing stuff")

# input for different tests
create_many_big_matrices <- function(nsubs=10, xdim=20, ydim=10, ...) {
    lapply(1:nsubs, function(x) as.big.matrix(scale(matrix(rnorm(xdim*ydim), xdim, ydim)), ...))
}

test_that("parts of compute_subdist_worker2 works", {
    nvoxs <- 100; nsubs <- 10
    funclist <- create_many_big_matrices(nsubs=nsubs, ydim=nvoxs)
    dmat <- big.matrix(nsubs^2, nvoxs, type="double", shared=FALSE)
    
    incols <- c(1,10)
    cormaps <- vbca_batch2(funclist, incols, ztransform=TRUE, shared=FALSE)
    
    i <- 2
    inds <- incols[1]:incols[2]
    voxs <- 1:nvoxs
    subsMap <- big.matrix(nvoxs-1, nsubs, type="double", shared=FALSE)
    
    # Combining and scaling
    .Call("CombineSubMapsMain", cormaps, subsMap@address, as.double(i), 
          as.double(voxs[-inds[i]]), as.double(nvoxs-1), as.double(nsubs))
    ref <- scale(sapply(cormaps, function(x) x[i,voxs[-inds[i]]]))
    expect_that(ref, is_equivalent_to(as.matrix(subsMap)))
    
    # Correlation
    big_cor(x=subsMap, z=dmat, z_firstCol=inds[i], z_lastCol=inds[i])
    m <- matrix(dmat[,inds[i]], nsubs, nsubs)
    mref <- cor(ref)
    expect_that(m, equals(mref))
    
    # Distance
    .Call("big_add_scalar", dmat, as.double(-1), as.double(1), 
            as.double(inds[i]), as.double(inds[i]));
    m <- matrix(dmat[,inds[i]], nsubs, nsubs)
    mref <- 1-cor(ref)
    expect_that(m, equals(mref))
})


library(connectir)
bm <- big.matrix(4, 4, init=2, type="double")
m <- .Call("test_func", bm)
m[,]
bm[,]



library(connectir)
create_many_big_matrices <- function(nsubs=10, xdim=20, ydim=10, ...) {
    lapply(1:nsubs, function(x) as.big.matrix(scale(matrix(rnorm(xdim*ydim), xdim, ydim)), ...))
}
voxinds <- 1:1000; nvoxs <- length(voxinds)
nsubs <- 20
funcs <- create_many_big_matrices(nsubs, 150, nvoxs, type="double", shared=FALSE)

# File-Backed
start.time <- Sys.time()
blocks <- niftir.split.indices(1, nvoxs, by=100)
dmatsA <- big.matrix(nsubs^2, nvoxs, type="double", backingpath=getwd(), 
                    backingfile="dmatsA.bin", descriptorfile="dmatsA.desc")
for (i in 1:blocks$n) {
    vcat(TRUE, "block %i", i)
    firstSeed <- blocks$start[i]; lastSeed <- blocks$ends[i]
    firstDist <- 1; lastDist <- lastSeed - firstSeed + 1
    compute_subdist2(funcs, firstSeed, lastSeed, 
                    dmatsA, firstDist, lastDist, 
                    blocksize=10, ztransform=TRUE, method="pearson", 
                    type="double")
    flush(dmatsA)
    dmatsA <- free.memory(dmatsA, getwd())
}
end.time <- Sys.time()
vcat(TRUE, "Took %.2f minutes\n", as.numeric(end.time-start.time, units="mins"))

# Non-File Backed
start.time <- Sys.time()
blocks <- niftir.split.indices(1, nvoxs, by=100)
dmatsB <- big.matrix(nsubs^2, nvoxs, type="double", shared=FALSE)
for (i in 1:blocks$n) {
    vcat(TRUE, "block %i", i)
    firstSeed <- blocks$start[i]; lastSeed <- blocks$ends[i]
    firstDist <- 1; lastDist <- lastSeed - firstSeed + 1
    compute_subdist2(funcs, firstSeed, lastSeed, 
                    dmatsB, firstDist, lastDist, 
                    blocksize=10, ztransform=TRUE, method="pearson", 
                    type="double")
}
tmp <- deepcopy(dmatsB, type="double", backingpath=getwd(), 
                    backingfile="dmatsB.bin", descriptorfile="dmatsB.desc")
end.time <- Sys.time()
vcat(TRUE, "Took %.2f minutes\n", as.numeric(end.time-start.time, units="mins"))







opts <- list(superblocksize=0, blocksize=0, memlimit=6, verbose=TRUE)
get_mdmr_memlimit(opts, nsubs=300, nvoxs=30000, nperms=15000, nfactors=4)

opts <- list(superblocksize=0, blocksize=0, memlimit=20, verbose=TRUE)
get_subdist_memlimit(opts, nsubs=300, nvoxs=30000, subs.ntpts=rep(300, 200))
