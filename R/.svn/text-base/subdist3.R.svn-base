.subdist_distance <- function(seedMaps, dmats, colind, 
                              transpose=FALSE, method="pearson", ...)
{
    if (method == "pearson") {
        scale_fast(seedMaps, to.copy=FALSE, byrows=transpose)
        .Call("subdist_pearson_distance", seedMaps, dmats, as.double(colind), 
              as.logical(transpose), PACKAGE="connectir")
    } else if (method == "shrink.pearson") {
        library(corpcor)
        if (transpose) seedMaps <- t(seedMaps[,])
        dmats[,colind] <- 1 - cor.shrink(seedMaps, verbose=FALSE, ...)[,]
    } else if (method == "icov") {
        library(glasso)
		# since we are just centering
		center_fast(seedMaps, to.copy=FALSE, byrows=transpose)
		# big_cor will give covariance & not correlation matrix
		oc <- big_cor(seedMaps, byrows=transpose)
		if (transpose)
			r <- norm_glasso(t(oc[,]), ...)
		else
			r <- norm_glasso(oc[,], ...)
        dmats[,colind] <- as.vector(1 - r)
		rm(oc, r); gc(F,T)
    } else {
        vstop("Unrecognized method %s", method)
    }
}

test_sdist <- function(...) .subdist_distance(...)

compute_subdist_wrapper3 <- function(sub.funcs1, list.dists, 
                                    blocksize, superblocksize, 
                                    sub.funcs2=sub.funcs2, 
                                    design_mat=NULL, 
                                    verbose=1, parallel=FALSE, ...)
{
    verbosity <- verbose
    verbose <- as.logical(verbose)
    sdist <- list.dists$sdist
    gdist <- list.dists$gdist
    bpath <- list.dists$bpath
    zcheck1 <- c(); zcheck2 <- c()
    
    nsubs <- length(sub.funcs1)
    if (nsubs != length(sub.funcs2))
        stop("length mismatch between 2 set of functional files")
    nvoxs <- ncol(sdist)
    superblocks <- niftir.split.indices(1, nvoxs, by=superblocksize)
    
    if (!is.null(design_mat)) {
        k <- qlm_rank(design_mat)
        if (k < ncol(design_mat))
            stop("design matrix is rank deficient")
    }
    
    vcat(verbose, "will run through %i large blocks", superblocks$n)
    for (i in 1:superblocks$n) {
        vcat(verbose, "large block %i", i)
        start.time <- Sys.time()
        
        firstSeed <- superblocks$start[i]; lastSeed <- superblocks$ends[i]
        firstDist <- 1; lastDist <- lastSeed - firstSeed + 1
        ncol <- lastDist
        
        # create temporary RAM-based matrix
        vcat(verbose, "...creating temporary distance matrices")
        tmp_sdist <- big.matrix(nsubs^2, ncol, type="double", shared=parallel)
        
        # subdist
        vcat(verbose, "...compute distances")
        compute_subdist3(sub.funcs1, firstSeed, lastSeed, sub.funcs2, 
                         tmp_sdist, firstDist, lastDist, 
                         blocksize=blocksize, design_mat=design_mat, 
                         verbose=verbosity, parallel=parallel, type="double", 
                         ...)
        
        # save subdist
        vcat(verbose, "...saving distances")
        firstCol <- superblocks$start[i]; lastCol <- superblocks$ends[i]
        sub_sdist <- sub.big.matrix(sdist, firstCol=firstCol, lastCol=lastCol, 
                                    backingpath=bpath)
        deepcopy(x=tmp_sdist, y=sub_sdist)
        ## checks
        tmp <- (sub_sdist[2,]!=0)*1 + 1
        zcheck1 <- c(zcheck1, tmp)
        if (any(tmp==1))
            vcat(verbose, "...there are %i bad voxels", sum(tmp==1))
        ## clear file-backed RAM usage
        flush(sub_sdist); flush(sdist)
        rm(sub_sdist); gc(FALSE, TRUE)
        sdist <- free.memory(sdist, bpath)

        # gower centered matrices
        vcat(verbose, "...gower centering")
        sub_gdist <- sub.big.matrix(gdist, firstCol=firstCol, lastCol=lastCol, 
                                    backingpath=bpath)
        gower.subdist2(tmp_sdist, outmat=sub_gdist, verbose=verbosity, parallel=parallel)
        ## checks
        tmp <- (sub_gdist[2,]!=0)*1 + 1
        zcheck2 <- c(zcheck2, tmp)
        if (any(tmp==1))
            vcat(verbose, "...there are %i bad voxels", sum(tmp==1))
        ## clear file-backed RAM usage
        flush(sub_gdist); flush(gdist)
        rm(sub_gdist); gc(FALSE, TRUE)
        gdist <- free.memory(gdist, bpath)
        
        # remove temporary matrix
        rm(tmp_sdist); gc(FALSE, TRUE)
        
        # how long?
        end.time <- Sys.time()
        time.total <- as.numeric(end.time-start.time, units="mins")
        time.left <- time.total*(superblocks$n-i)
        vcat(verbose, "...took %.1f minutes (%.1f minutes left)\n", 
             time.total, time.left)
    }
    
    list(sdist=zcheck1, gdist=zcheck2)
}

compute_subdist3 <- function(sub.funcs1, firstSeed, lastSeed, sub.funcs2, 
                             dmats, firstDist, lastDist, 
                             blocksize=floor(ncol(dmats)/getDoParWorkers()), 
                             design_mat=NULL, verbose=1, parallel=FALSE, 
                             ...)
{
    nseeds <- lastSeed - firstSeed + 1
    ndists <- lastDist - firstDist + 1
    if (nseeds != ndists)
        stop("length mismatch between # of seeds  and # of distance matrices")
    seeds <- firstSeed:lastSeed
    dists <- firstDist:lastDist
    
    blocks <- niftir.split.indices(2, nseeds-1, by=blocksize)
    use_shared <- ifelse(parallel, TRUE, FALSE)
    progress <- ifelse(as.logical(verbose), "text", "none")
    inform <- verbose==2
    verbose <- as.logical(verbose)
    
    if (!is.big.matrix(sub.funcs1[[1]]) || !is.big.matrix(sub.funcs2[[1]]) || !is.big.matrix(dmats))
        stop("inputs and outputs must be big matrices")
    if (parallel && (!is.shared(sub.funcs1[[1]]) || !is.shared(sub.funcs2[[1]]) || !is.shared(dmats)))
        stop("if running in parallel inputs and outputs must be of type shared")
    
    if (is.null(design_mat)) {
        dfun <- function(starti, lasti, ...) {
            sub.firstSeed <- seeds[starti]; sub.lastSeed <- seeds[lasti]
            sub.firstDist <- dists[starti]; sub.lastDist <- dists[lasti]
            compute_subdist_worker3(sub.funcs1, sub.firstSeed, sub.lastSeed, sub.funcs2, 
                                    dmats, sub.firstDist, sub.lastDist, 
                                    shared=FALSE, ...)
            return(NULL)
        }
    } else {
        dfun <- function(starti, lasti, ...) {
            sub.firstSeed <- seeds[starti]; sub.lastSeed <- seeds[lasti]
            sub.firstDist <- dists[starti]; sub.lastDist <- dists[lasti]
            compute_subdist_worker3_regress(sub.funcs1, sub.firstSeed, sub.lastSeed, sub.funcs2, 
                                            dmats, sub.firstDist, sub.lastDist, 
                                            design_mat, 
                                            shared=FALSE, ...)
            return(NULL)
        }
    }
    
    # Test
    vcat(verbose, "...running a test on first seed")
    dfun(1, 1, ...)
    check_dmat(matrix(dmats[,dists[1]], sqrt(nrow(dmats))))
    vcat(verbose, "...running a test on last seed")
    dfun(ndists, ndists, ...)
    check_dmat(matrix(dmats[,dists[ndists]], sqrt(nrow(dmats))))
    
    # Subdist Calculation
    vcat(verbose, "...now the real deal with %i blocks and %i seeds", blocks$n, nseeds-2)
    llply(1:blocks$n, function(i, ...) {
        starti <- blocks$starts[i]; lasti <- blocks$ends[i]
        dfun(starti, lasti, ...)
    }, ..., .progress=progress, .parallel=parallel, .inform=inform)
    
    invisible(TRUE)
}

compute_subdist_worker3 <- function(sub.funcs1, firstSeed, lastSeed, sub.funcs2, 
                                    dmats, firstDist, lastDist, 
                                    ztransform=FALSE, method="pearson", 
                                    type="double", shared=FALSE, ...)
{
    nsubs <- length(sub.funcs1)
    nvoxs <- ncol(sub.funcs2[[1]])
    nseeds <- lastSeed - firstSeed + 1
    ndists <- lastDist - firstDist + 1
    if (nseeds != ndists)
        stop("mismatch in length of seed and distance matrix indices")
    voxs <- 1:nvoxs
    seeds <- firstSeed:lastSeed
    dists <- firstDist:lastDist
    
    subs.cormaps <- vbca_batch3(sub.funcs1, c(firstSeed, lastSeed), sub.funcs2, 
                                ztransform=ztransform, 
                                type=type, shared=shared, ...)
    
    # check first subject for any inf or NAs and save those
    
    seedCorMaps <- big.matrix(nvoxs, nsubs, type=type, shared=shared)
    for (i in 1:nseeds) {
        .Call("subdist_combine_submaps", subs.cormaps, as.double(i), 
              as.double(voxs), seedCorMaps, PACKAGE="connectir")
        .subdist_distance(seedCorMaps, dmats, dists[i], FALSE, method)
    }
    
    rm(subs.cormaps, seedCorMaps)
    gc(FALSE, TRUE)
    
    return(dmats)
}

compute_subdist_worker3_regress <- function(sub.funcs1, firstSeed, lastSeed, sub.funcs2, 
                                            dmats, firstDist, lastDist, 
                                            design_mat, 
                                            ztransform=FALSE, method="pearson", 
                                            type="double", shared=FALSE, ...)
{
    nsubs <- length(sub.funcs1)
    nvoxs <- ncol(sub.funcs2[[1]])
    nseeds <- lastSeed - firstSeed + 1
    ndists <- lastDist - firstDist + 1
    if (nseeds != ndists)
        stop("mismatch in length of seed and distance matrix indices")
    voxs <- 1:nvoxs
    seeds <- firstSeed:lastSeed
    dists <- firstDist:lastDist
    
    subs.cormaps <- vbca_batch3(sub.funcs1, c(firstSeed, lastSeed), sub.funcs2, 
                                ztransform=ztransform, 
                                type=type, shared=shared, ...)
    
    seedCorMaps <- big.matrix(nsubs, nvoxs, type=type, shared=shared)
    r_seedCorMaps <- big.matrix(nsubs, nvoxs, type=type, shared=shared)
    for (i in 1:nseeds) {
        .Call("subdist_combine_and_trans_submaps", subs.cormaps, as.double(i), 
              as.double(voxs), seedCorMaps, PACKAGE="connectir")
        qlm_residuals(seedCorMaps, design_mat, FALSE, r_seedCorMaps, TRUE)
        .subdist_distance(r_seedCorMaps, dmats, dists[i], TRUE, method)
    }
    
    rm(subs.cormaps, seedCorMaps)
    gc(FALSE, TRUE)
    
    return(dmats)
}
