# threshType: 0=mean all; 1=thresh&mean; 2=thresh&bin&sum...
gcor <- function(bigmat, blocksize, ztransform=FALSE, 
                thresh=0, threshType=0, verbosity=1, 
                parallel=FALSE, ...) 
{
    ## input info
    nvoxs <- ncol(bigmat)
    vox_inds <- as.double(1:nvoxs)
    if (ztransform)
        thresh <- atanh(thresh)
    blocks <- niftir.split.indices(1, nvoxs, by=blocksize)
    thresh <- as.double(thresh)
    threshType <- as.integer(threshType)
    
    verbose <- as.logical(verbosity)
    inform <- verbosity==2
    progress <- ifelse(verbose, "text", "none")
    
    ## function that computes correlation maps for subset of voxels
    ## + gets average of each correlation map
    dfun <- function(i) {
        inds_CHUNK <- vox_inds[c(blocks$starts[i],blocks$ends[i])]
        cormat_CHUNK <- vbca2(bigmat, inds_CHUNK, ztransform=ztransform, ...)
        
        inds <- as.double(inds_CHUNK[1]:inds_CHUNK[2])
        gs <- .Call("gcor_worker", cormat_CHUNK, thresh, threshType, 
                    inds, vox_inds, package="connectir")
        rm(inds_CHUNK, cormat_CHUNK); gc(FALSE, TRUE)
        
        as.vector(gs)
    }
    
    gcor <- llply(1:blocks$n, dfun, .progress=progress, .inform=inform, 
                  .parallel=parallel)
    gcor <- unlist(gcor)
    if (ztransform && (threshType == 0 || threshType == 1 || threshType == 3))
        gcor <- tanh(gcor)
    
    return(gcor)
}

# This is a simple version of kendall
kendall_ref <- function(ratings) {
    ratings <- as.matrix(na.omit(ratings))
    ns <- nrow(ratings)
    nr <- ncol(ratings)
    
    ratings.rank <- apply(ratings, 2, rank)
    coeff <- (12 * var(rowSums(ratings.rank)) * (ns - 
        1))/(nr^2 * (ns^3 - ns))
    
    return(coeff)
}

# This computes a kendall's W examining the consistency of
# each voxel's connectivity map across participants
kendall <- function(subs.bigmats, blocksize, ztransform=FALSE, parallel=FALSE, 
                    verbosity=1) 
{
    nsubs <- length(subs.bigmats)
    nvoxs <- ncol(subs.bigmats[[1]])
    voxs <- as.double(1:nvoxs)
    blocks <- niftir.split.indices(1, nvoxs, by=blocksize)
    inform <- verbosity==2
    progress <- ifelse(verbosity>0, "text", "none")
    
    kfun <- function(i) {
        incols <- voxs[c(blocks$starts[i],blocks$ends[i])]
        cormats <- vbca_batch2(subs.bigmats, incols, ztransform=ztransform, 
                               type="double", shared=FALSE)
        seeds <- as.double(incols[1]:incols[2])
        seedMaps <- big.matrix(nvoxs-1, nsubs, type="double", shared=FALSE)       
        coeffs <- .Call("voxelwise_kendall", cormats, seedMaps, seeds, voxs)
        rm(cormats, seedMaps); gc(FALSE, TRUE)
        return(coeffs)
    }
    
    ws <- llply(1:blocks$n, kfun, .progress=progress, .parallel=parallel, .inform=inform)
    ws <- unlist(ws)
    
    return(ws)
}

kendall3 <- function(subs.func1, blocksize, subs.func2=subs.func1, 
                     design_mat=NULL, ztransform=FALSE, 
                     parallel=FALSE, verbose=TRUE) 
{
    nsubs <- length(subs.func1)
    if (is.null(subs.func2)) subs.func2 <- subs.func1
    if (nsubs != length(subs.func2))
        stop("length mismatch in 2 set of functionals")
    nseeds <- ncol(subs.func1[[1]])
    seeds <- as.double(1:nseeds)
    blocks <- niftir.split.indices(1, nseeds, by=blocksize)
    nvoxs <- ncol(subs.func2[[1]])
    voxs <- as.double(1:nvoxs)
    if (!is.null(design_mat)) {
        if (!is.big.matrix(design_mat))
            stop("design_mat must be a big matrix")
        if (nrow(design_mat) != nsubs)
            stop("design mat has row mismatch (wrong # of subjects)")
    }
    inform <- verbose
    progress <- ifelse(verbose, "text", "none")

    if (is.null(design_mat)) {
        kfun <- function(i) {
            incols <- seeds[c(blocks$starts[i],blocks$ends[i])]
            cormats <- vbca_batch3(subs.func1, incols, subs.func2, ztransform=ztransform, 
                                   type="double", shared=FALSE)
            seeds <- as.double(incols[1]:incols[2])
            seedMaps <- big.matrix(nvoxs, nsubs, type="double", shared=FALSE)       
            coeffs <- .Call("voxelwise_kendall3", cormats, seedMaps, seeds, voxs)
            rm(cormats, seedMaps); gc(FALSE, TRUE)
            return(coeffs)
        }
    } else {
        kfun <- function(i) {
            incols <- seeds[c(blocks$starts[i],blocks$ends[i])]
            cormats <- vbca_batch3(subs.func1, incols, subs.func2, ztransform=ztransform, 
                                   type="double", shared=FALSE)
            seeds <- as.double(incols[1]:incols[2])
            seedMaps <- big.matrix(nvoxs, nsubs, type="double", shared=FALSE)       
            coeffs <- .Call("voxelwise_kendall3_regress", cormats, seedMaps, 
                            design_mat, seeds, voxs)
            rm(cormats, seedMaps); gc(FALSE, TRUE)
            return(coeffs)
        }
    }
        
    ws <- llply(1:blocks$n, kfun, .progress=progress, .parallel=parallel, .inform=inform)
    ws <- unlist(ws)
    
    return(ws)
}

reho_worker <- function(mat, ...) {
    return(.Call("kendall_worker", mat, ...))
}

reho <- function(bigmat, nei=1, nei.dist=3, min.nei=0.5, verbose=TRUE, 
                 parallel=FALSE, FUN=reho_worker, ...)
{
    header <- bigmat@header
    header$dim <- header$dim[1:3]
    header$pixdim <- header$pixdim[1:3]
    progress <- ifelse(verbose, "text", "none")
    
    mask <- bigmat@mask
    nvoxs <- ncol(bigmat)
    if (sum(mask) != nvoxs)
        stop("Number of TRUE elems in mask does not equal number of columns in bigmat")
    
    dims <- header$dim
    moffsets <- expand.grid(list(i=-nei:nei, j=-nei:nei, k=-nei:nei))
    dist <- rowSums(abs(moffsets))
    moffsets <- moffsets[dist<=nei.dist,]
    offsets <- moffsets$k*dims[1]*dims[2] + moffsets$j*dims[1] + moffsets$i
    
    mat2arr_inds <- which(mask)
    arr2mat_inds <- mask*1
    arr2mat_inds[mask] <- 1:nvoxs
    
    opts <- list()
    opts$min.nei <- ceiling(nrow(moffsets)*min.nei)
    opts$rawi.min <- 0
    opts$rawi.max <- length(mask)+1
    opts$dim <- header$dim
    opts$mat2arr_inds <- mat2arr_inds
    opts$arr2mat_inds <- arr2mat_inds
    opts$offsets <- offsets
    
    rfun <- function(i, ...) {
        raw_inds <- opts$offsets + opts$mat2arr_inds[i]
        raw_inds <- raw_inds[raw_inds > opts$rawi.min & raw_inds < opts$rawi.max]
        filt_inds <- opts$arr2mat_inds[raw_inds]
        filt_inds <- filt_inds[filt_inds>0]
        
        if (length(filt_inds) < opts$min.nei)
            return(0)
        
        x <- deepcopy(bigmat, cols=filt_inds)
        rs <- FUN(x, ...)
        rm(x); gc(FALSE, TRUE)
        
        return(rs)
    }
    
    reho.vals <- llply(1:nvoxs, rfun, ..., 
                       .progress=progress, .inform=TRUE, .verbose=verbose)
    reho.vals <- unlist(reho.vals)
    
    return(list(reho=reho.vals, hdr=header, mask=mask))
}

