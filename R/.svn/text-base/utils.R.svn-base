vcat <- function(verbose, msg, ..., newline=TRUE) {
    if (verbose) {
        cat(sprintf(msg, ...))
        if (newline) cat("\n")
    }
}

vstop <- function(msg, ...) stop(sprintf(msg, ...))

vsystem <- function(verbose, cmd, ...) {
    vcat(verbose, cmd, ...)
    ret <- system(sprintf(cmd, ...))
    if (ret != 0)
        stop("command failed")
}

set_parallel_procs <- function(nforks=1, nthreads=1, verbose=FALSE, force=FALSE) {
    vcat(verbose, "Setting %i parallel forks", nforks)
    suppressPackageStartupMessages(library("doMC"))
    registerDoMC()
    nprocs <- getDoParWorkers()
    if (nforks > nprocs) {
        msg <- sprintf("# of forks %i is greater than the actual # of processors (%i)", 
                          nforks, nprocs)
        if (force == TRUE) {
            warning(msg, immediate.=TRUE)
        } else {
            vstop(msg)
        }
    }
    options(cores=nforks)
    
    vcat(verbose, "Setting %i threads for matrix algebra operations", 
         nthreads)
    #nprocs <- omp_get_max_threads()
    if (nthreads > nprocs) {
        msg <- sprintf("# of threads %i is greater than the actual # of processors (%i)", 
                          nthreads, nprocs)
        if (force == TRUE) {
          warning(msg, immediate.=TRUE)
        } else {
          vstop(msg)
        }
    }
    
    if (existsFunction("setMKLthreads", where=topenv(.GlobalEnv))) {
        vcat(verbose, "...using Intel's MKL")
        setMKLthreads(nthreads)
    } else {
        # cover all our blases
        vcat(verbose, "...using GOTOBLAS or Other")
        suppressPackageStartupMessages(library("blasctl"))
        blas_set_num_threads(nthreads)
        omp_set_num_threads(nthreads)
    }
    
    # Not sure if these env vars matter?
    Sys.setenv(OMP_NUM_THREADS=nthreads)
    Sys.setenv(GOTO_NUM_THREADS=nthreads)
    Sys.setenv(MKL_NUM_THREADS=nthreads)
    Sys.setenv(OMP_DYNAMIC=TRUE)
    Sys.setenv(MKL_DYNAMIC=TRUE)
    
    invisible(TRUE)
}

# dims and mask useful is pass vector that represents 3D dataset
cluster.table <- function(x, vox.thr=0, dims=NULL, mask=NULL, nei=1, nei.dist=3, pad=1) {
    if (is.null(mask))
        mask <- rep(T, length(x))
    if (is.null(dims)) {
        if (length(dim(x)) != 3)
            stop("If dims isn't provided, then x must be a 3D array")
        dims <- dim(x)
    }
    
    nx <- dims[1]; ny <- dims[2]; nz <- dims[3]
    if (pad > 0) {
        mask.pad <- array(F, c(nx+pad*2, ny+pad*2, nz+pad*2))
        mask.pad[(pad+1):(pad+nx), (pad+1):(pad+ny), (pad+1):(pad+nz)] <- T
        nx <- dim(mask.pad)[1]; ny <- dim(mask.pad)[2]; nz <- dim(mask.pad)[3]
        mask.pad <- as.vector(mask.pad)
        tmp <- rep(F, length(mask.pad))
        tmp[mask.pad] <- mask
        mask <- tmp
        rm(tmp)
    }
    
    xfull <- vector("numeric", nx*ny*nz)
    mfull <- vector("logical", nx*ny*nz)
    
    # coords <- expand.grid(list(x=1:nx, y=1:ny, z=1:nz))
    
    # Get neighbours to check for clusters
    # default is 27 (face, edge, corner touching)
    nmat <- expand.grid(list(i=-nei:nei, j=-nei:nei, k=-nei:nei))
    voxdists <- rowSums(abs(nmat))
    nmat <- nmat[voxdists<=nei.dist,]
    offsets <- nmat$k*nx*ny + nmat$j*nx + nmat$i
    rm(nmat)
    
    # threshold
    tmp.mask <- x > vox.thr
    nvoxs <- sum(tmp.mask)
    mfull[mask] <- tmp.mask
    xfull[mask][tmp.mask] <- x[tmp.mask]
    rm(tmp.mask)
    
    # Create a list of indices
    mm <- rep(T, nvoxs)
    mat2arr_inds <- which(mfull)
    arr2mat_inds <- mfull*1
    arr2mat_inds[mfull] <- 1:nvoxs
    maxi <- length(mfull) + 1
    
    cfull <- mfull*0
    m <- mfull
    sizes <- c()
    masses <- c()
    n <- nvoxs
    nc <- 0
    while (n > 0) {
        i <- mat2arr_inds[mm][1]
        
        clust.n <- 1
        clust.mass <- xfull[i]
        n <- n - 1
        m[i] <- F
        mm[mm][1] <- F
        nc <- nc + 1
        cfull[i] <- nc
        
        inds <- offsets + i
        inds <- inds[inds > 0 & inds < maxi]
        tmp_inds <- inds[m[inds]]
        
        if (length(tmp_inds) > 0) {
            m[tmp_inds] <- F
            mm[arr2mat_inds[tmp_inds]] <- F

            while(length(tmp_inds) > 0) {
                i <- tmp_inds[1]
                inds <- offsets + i
                inds <- inds[inds > 0 & inds < maxi]
                inds <- inds[m[inds]]
                if (length(inds) > 0) {
                    m[inds] <- F
                    mm[arr2mat_inds[inds]] <- F
                }
                tmp_inds <- c(tmp_inds[-1], inds)
                
                clust.n <- clust.n + 1
                clust.mass <- clust.mass + xfull[i]
                n <- n - 1
                cfull[i] <- nc
            }
        }
        
        masses <- c(clust.mass, masses)
        sizes <- c(clust.n, sizes)
    }
    
    if (pad > 0)
        cfull <- cfull[mask.pad]
    dim(cfull) <- dims
    
    if (length(sizes) == 0) {
        nc <- 0
        sizes <- 0
        masses <- 0
    }
    
    # RETURNS:
    # max.size: maximum cluster size
    # max.mass: maximum cluster mass
    # size: different cluster sizes
    # mass: different cluster masses
    # clust: 3D array like x but with clusters numbered
    #        (cluster # here related to index of size and mass)
    list(
        nclust=nc, 
        max.size=max(sizes), 
        max.mass=max(masses), 
        max.rel=max(masses/sizes), 
        size=sizes,
        mass=masses,
        rel=masses/sizes, 
        clust=cfull
    )
}

big_rowsum <- function(bigmat, firstCol=1, lastCol=ncol(bigmat)) {
    as.vector(.Call("big_rowsum", bigmat, as.double(firstCol), as.double(lastCol), 
              PACKAGE="connectir"))
}

big_rowmean <- function(bigmat, firstCol=1, lastCol=ncol(bigmat)) {
    as.vector(.Call("big_rowmean", bigmat, as.double(firstCol), as.double(lastCol), 
              PACKAGE="connectir"))
}

# x = a*x + b
big_add_multiply_scalar <- function(x, a=1, b=0, firstCol=1, lastCol=ncol(x)) {
    .Call("big_add_multiply_scalar", x, as.double(a), as.double(b), 
          as.double(firstCol), as.double(lastCol), PACKAGE="connectir")
    invisible(x)
}
