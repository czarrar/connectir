check_dmat <- function(dmat) {
    TOL <- .Machine$double.eps ^ 0.5
    dmat <- abs(as.matrix(dmat))
    diag_dmat <- diag(dmat)
    off_dmat <- dmat[lower.tri(dmat)]
    
    if (any(is.na(dmat)))
        stop("NAs were present in distance matrix")
    
    if (all(dmat < TOL))
        stop("All zeros in distance matrix")
    
    if (all(diag_dmat>TOL)) {
        warning("Diagonal of distance matrix is all non-zeros\n")
    } else if (any(diag_dmat>TOL)) {
        cat("Diagonal of distance matrix has non-zeros\n")
        print(diag_dmat)
    }
    
    if (any(off_dmat<TOL))
        cat("Off-diagonal of distance matrix has some zeros\n")
}

check_gmat <- function(gmat) {
    TOL <- .Machine$double.eps ^ 0.5
    dmat <- abs(as.matrix(gmat))
    
    if (any(is.na(gmat)))
        stop("NAs were present in gower's centered distance matrix")
    if (all(dmat < TOL))
        stop("All zeros in gower's centered distance matrix")
}

# This checks that everything in folder is good
check_subdist <- function(sdir) {
    # Checking Functions
    checkpath <- function(p, is.dir) { 
        if (!file.exists(p)) {
            if (is.dir)
                stop("Directory: ", p, " does not exist but required")
            else
                stop("File: ", p, " does not exist but required")
        }
    }
    checkthing <- function(comparison, ...) {
        if (!comparison)
            stop(...)
    }
    
    # File/directory paths to check
    sdir <- abspath(sdir)
    optsfile <- file.path(sdir, "options.rda")
    infuncdir <- file.path(sdir, "input_funcs")
    maskfile <- file.path(sdir, "mask.nii.gz")
    distfiles <- file.path(sdir, 
        c("subdist.desc", "subdist.bin", "subdist_gower.desc", "subdist_gower.bin")
    )
    sdistfile <- file.path(sdir, "subdist.desc")
    gdistfile <- file.path(sdir, "subdist_gower.desc")
    
    # Check main directory
    checkpath(sdir, TRUE)
    
    # Check opts and read them in
    checkpath(optsfile, FALSE)
    opts <- NULL
    load(optsfile)
    if (is.null(opts))
        stop("Optsfile doesn't have opts variable!")
    
    # Check input funcs
    checkpath(infuncdir, TRUE)
    checkthing(
        length(list.files(infuncdir))==length(opts$infiles), 
        "Missing some input functional files"
    )
    
    # Check input masks
    checkpath(maskfile, FALSE)
    
    # Check subdist and related files
    lapply(distfiles, checkpath, FALSE)
    
    # Count number of seed voxels
    mask <- read.mask(maskfile)
    nvoxs <- sum(mask)
    
    # Read in subdist and seedmask + check
    lapply(c(sdistfile, gdistfile), function(f) {
        x <- attach.big.matrix(f)
        ## make sure appropriate number of voxels
        checkthing(
            nvoxs==ncol(x),
            "# of seed voxels does not match # of columns in ", f
        )
        ## mask sure appropriate number of subjects
        nsubs <- sqrt(nrow(x))
        checkthing(
            length(opts$infiles)==nsubs,
            "# of subjects does not match # of rows in ", f
        )
    })
}

# This creates a new subdist directory and relevant files
# and returns a new subdist object
create_subdist <- function(outdir, infiles1, mask1, infiles2, mask2, opts, ...) {
    if (file.exists(outdir))
        stop("output directory cannot exist")
    
    if (is.null(infiles2) || is.null(mask2))
        use.set2 <- FALSE
    else
        use.set2 <- TRUE
    
    infuncdir <- file.path(outdir, "input_funcs")
    
    if (use.set2)
        infuncdir2 <- file.path(outdir, "input_funcs2")
    
    # Create directories
    vcat(opts$verbose, "...creating directories in %s", outdir)
    dir.create(outdir)
    dir.create(infuncdir)
    if (use.set2)
        dir.create(infuncdir2)
    
    # Create symlinks for the input funcs
    if (!opts$"no-link-functionals") {
        vcat(opts$verbose, "...creating soft links to subject functional data")
        for (i in 1:length(infiles1)) {
            from <- infiles1[i]
            to <- file.path(infuncdir, sprintf("scan%04i.%s", i, getext(from)))
            file.symlink(from, to)
        }
        if (use.set2) {
            for (i in 1:length(infiles2)) {
                from <- infiles2[i]
                to <- file.path(infuncdir2, sprintf("scan%04i.%s", i, getext(from)))
                file.symlink(from, to)
            }
        }
    }
    
    # Get a header file from the first functional
    hdr <- read.nifti.header(infiles1[1])
    if (length(hdr$dim) == 4) {
        hdr$dim <- hdr$dim[1:3]
        hdr$pixdim <- hdr$pixdim[1:3]
    } else if (length(hdr$dim) == 2) {
        hdr$dim <- c(hdr$dim[2], 1)
        hdr$pixdim <- c(hdr$pixdim[2], 1)
    }
    
    if (use.set2) {
        hdr2 <- read.nifti.header(infiles2[1])
        if (length(hdr2$dim) == 4) {
            hdr2$dim <- hdr2$dim[1:3]
            hdr2$pixdim <- hdr2$pixdim[1:3]
        } else if (length(hdr2$dim) == 2) {
            hdr2$dim <- c(hdr2$dim[2], 1)
            hdr2$pixdim <- c(hdr2$pixdim[2], 1)
        }
    }
    
    # Write the brain masks
    vcat(opts$verbose, "...saving masks")
    outfile <- file.path(outdir, "mask.nii.gz")
    write.nifti(mask1, hdr, outfile=outfile, odt="char")
    if (use.set2) {
        outfile <- file.path(outdir, "mask2.nii.gz")
        write.nifti(mask2, hdr2, outfile=outfile, odt="char")
    }
    
    # Copy over standard brain
    vcat(opts$verbose, "...copying background image")
    if (!is.null(opts$bg))
        file.copy(opts$bg, file.path(outdir, "bg_image.nii.gz"))
    
    # Save options
    vcat(opts$verbose, "...saving options")
    opts$outdir <- outdir
    opts$infiles <- infiles1
    opts$infiles2 <- infiles2
    save(opts, file=file.path(outdir, "options.rda"))
    
    # Create file-backed subject distances and gower matrices
    vcat(opts$verbose, "...creating file-backed distance matrices")
    nsubs <- length(infiles1)
    nvoxs <- sum(mask1)
    sdist <- big.matrix(nsubs^2, nvoxs, type="double", 
                        backingpath=outdir, 
                        backingfile="subdist.bin", 
                        descriptorfile="subdist.desc")
    gdist <- big.matrix(nsubs^2, nvoxs, type="double", 
                        backingpath=outdir, 
                        backingfile="subdist_gower.bin", 
                        descriptorfile="subdist_gower.desc")
    
    # Create temporary subject distances matrix
    
    list(sdist=sdist, gdist=gdist, bpath=outdir)
}

compute_subdist_wrapper <- function(sub.funcs, list.dists, 
                                    blocksize, superblocksize, 
                                    design_mat=NULL, 
                                    verbose=1, parallel=FALSE, ...)
{
    verbosity <- verbose
    verbose <- as.logical(verbose)
    sdist <- list.dists$sdist
    gdist <- list.dists$gdist
    bpath <- list.dists$bpath
    zcheck1 <- c(); zcheck2 <- c()
    
    nsubs <- length(sub.funcs)
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
        compute_subdist2(sub.funcs, firstSeed, lastSeed, 
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

compute_subdist2 <- function(sub.funcs, firstSeed, lastSeed, 
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
    
    if (!is.big.matrix(sub.funcs[[1]]) || !is.big.matrix(dmats))
        stop("inputs and outputs must be big matrices")
    if (parallel && (!is.shared(sub.funcs[[1]]) || !is.shared(dmats)))
        stop("if running in parallel inputs and outputs must be of type shared")
    
    if (is.null(design_mat)) {
        dfun <- function(starti, lasti, ...) {
            sub.firstSeed <- seeds[starti]; sub.lastSeed <- seeds[lasti]
            sub.firstDist <- dists[starti]; sub.lastDist <- dists[lasti]
            compute_subdist_worker2(sub.funcs, sub.firstSeed, sub.lastSeed, 
                                    dmats, sub.firstDist, sub.lastDist, 
                                    shared=FALSE, ...)
            return(NULL)
        }
    } else {
        dfun <- function(starti, lasti, ...) {
            sub.firstSeed <- seeds[starti]; sub.lastSeed <- seeds[lasti]
            sub.firstDist <- dists[starti]; sub.lastDist <- dists[lasti]
            compute_subdist_worker2_regress(sub.funcs, sub.firstSeed, sub.lastSeed, 
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

compute_subdist_worker2 <- function(sub.funcs, firstSeed, lastSeed, 
                                    dmats, firstDist, lastDist, 
                                    ztransform=FALSE, method="pearson", 
                                    type="double", shared=FALSE, ...)
{
    nsubs <- length(sub.funcs)
    nvoxs <- ncol(sub.funcs[[1]])
    nseeds <- lastSeed - firstSeed + 1
    ndists <- lastDist - firstDist + 1
    if (nseeds != ndists)
        stop("mismatch in length of seed and distance matrix indices")
    voxs <- 1:nvoxs
    seeds <- firstSeed:lastSeed
    dists <- firstDist:lastDist
    
    subs.cormaps <- vbca_batch2(sub.funcs, c(firstSeed, lastSeed), 
                                ztransform=ztransform, 
                                type=type, shared=shared, ...)
    
    seedCorMaps <- big.matrix(nvoxs-1, nsubs, type=type, shared=shared, ...)
    for (i in 1:nseeds) {
        .Call("subdist_combine_submaps", subs.cormaps, as.double(i), 
              as.double(voxs[-seeds[i]]), seedCorMaps, PACKAGE="connectir")
        .subdist_distance(seedCorMaps, dmats, dists[i], FALSE, method)
    }
    
    rm(subs.cormaps, seedCorMaps)
    gc(FALSE, TRUE)
    
    return(dmats)
}

compute_subdist_worker2_regress <- function(sub.funcs, firstSeed, lastSeed, 
                                            dmats, firstDist, lastDist, 
                                            design_mat, 
                                            ztransform=FALSE, method="pearson", 
                                            type="double", shared=FALSE, ...)
{
    nsubs <- length(sub.funcs)
    nvoxs <- ncol(sub.funcs[[1]])
    nseeds <- lastSeed - firstSeed + 1
    ndists <- lastDist - firstDist + 1
    if (nseeds != ndists)
        stop("mismatch in length of seed and distance matrix indices")
    voxs <- 1:nvoxs
    seeds <- firstSeed:lastSeed
    dists <- firstDist:lastDist
    
    subs.cormaps <- vbca_batch2(sub.funcs, c(firstSeed, lastSeed), 
                                ztransform=ztransform, 
                                type=type, shared=shared, ...)
    
    seedCorMaps <- big.matrix(nsubs, nvoxs-1, type=type, shared=shared, ...)
    r_seedCorMaps <- big.matrix(nsubs, nvoxs-1, type=type, shared=shared, ...)
    for (i in 1:nseeds) {
        .Call("subdist_combine_and_trans_submaps", subs.cormaps, as.double(i), 
              as.double(voxs[-seeds[i]]), seedCorMaps, PACKAGE="connectir")
        qlm_residuals(seedCorMaps, design_mat, FALSE, r_seedCorMaps, TRUE)
        .subdist_distance(r_seedCorMaps, dmats, dists[i], TRUE, method)
    }
    
    rm(subs.cormaps, seedCorMaps)
    gc(FALSE, TRUE)
    
    return(dmats)
}

# bigmat: rows=subject distances, cols=voxels
# not that each column is a vectorized version of n x n matrix comparing subjects where n^2 = # of row elements
gower.subdist <- function(bigmat, gower.bigmat=NULL, verbose=TRUE, ...) {
    # this all does
    # G <- -0.5 * -(dmat*dmat) %*% (I - ones %*% t(ones)/n)
    
    # setup
    nc <- ncol(bigmat)
    n <- sqrt(nrow(bigmat))
    I <- diag(n)
    ones <- matrix(1, nrow=n)
    if (is.null(gower.bigmat))
        gower.bigmat <- deepcopy(bigmat, ...)
    
    # bigmat <- bigmat * bigmat
    pow(gower.bigmat, 2)
    
    # bigmat <- -(bigmat)/2
    #big_add_multiply_scalar(x=gower.bigmat, a=-0.5)
    dscal(ALPHA=-0.5, Y=gower.bigmat)
    
    # I - ones %*% t(ones)/n
    adj <- I - ones %*% t(ones)/n
    
    if (verbose)
        pb <- progressbar(nc)
    
    # newmat <- bigmat %*% adj
    gfun <- function(i) {
        if (verbose)
            update(pb, i)
        gower.vox <- matrix(gower.bigmat[,i], n, n)
        gower.bigmat[,i] <- as.vector(gower.vox %*% adj)
        return(NULL)
    }
    for (i in seq_len(nc))
        gfun(i)
    #TODO: seems like the parallel thing here takes longer than it should
    #if (getDoParWorkers() == 1) {
    #    for (i in seq_len(nc)) gfun(i)
    #} else {
    #    foreach(i = seq_len(nc)) %dopar% gfun(i)
    #}
    
    if (verbose)
        end(pb)
    
    return(gower.bigmat)
}

gower.subdist2 <- function(inmat, outmat=NULL, blocksize=ncol(inmat), 
                           verbose=1, parallel=FALSE, ...) 
{
    nr <- nrow(inmat)
    nc <- ncol(inmat)
    blocksize <- max(floor(blocksize/getDoParWorkers()), 1)
    blocks <- niftir.split.indices(1, nc, by=blocksize)
    use_shared <- ifelse(parallel, TRUE, FALSE)
    progress <- ifelse(as.logical(verbose), "text", "none")
    inform <- verbose==2
    verbose <- as.logical(verbose)
    
    if (is.null(outmat)) {
        outmat <- big.matrix(nr, nc, type="double", ...)
    } else if (nr != nrow(outmat) || nc != ncol(outmat)) {
        vstop("size mismatch between input (%ix%i) and output (%ix%i)", 
                nr, nc, nrow(outmat), ncol(outmat))
    }
    
    if (!is.big.matrix(inmat) || !is.big.matrix(outmat))
        stop("input and output must be big matrices")
    if (parallel && (!is.shared(inmat) || !is.shared(outmat)))
        stop("if running in parallel input and output must be of type shared")
    
    llply(1:blocks$n, function(i) {
        si <- blocks$starts[i]
        ei <- blocks$ends[i]
        .Call("big_gower", inmat, outmat, as.double(si), as.double(ei), 
              as.double(si), as.double(ei), PACKAGE="connectir")
    }, .progress=progress, .parallel=parallel, .inform=inform)
    
    return(outmat)
}

square.subdist <- function(bigmat, square.bigmat=NULL, ...) {
    # this all does
    # Amat <- -0.5 * -(dmat*dmat)
    
    # setup
    if (is.null(square.bigmat))
        square.bigmat <- deepcopy(bigmat, ...)
    
    # bigmat <- bigmat * bigmat
    .Call("BigPowMain", square.bigmat@address, as.double(2))
    
    # bigmat <- -(bigmat)/2
    dscal(ALPHA=-0.5, Y=square.bigmat)
    
    return(square.bigmat)
}

filter_subdist <- function(bigmat, subs=1:sqrt(nrow(bigmat)), voxs=1:ncol(bigmat), ...) {
    matinds <- matrix(1:nrow(bigmat), sqrt(nrow(bigmat)))
    matinds <- as.vector(matinds[subs,subs])
    deepcopy(bigmat, cols=voxs, rows=matinds, ...)
}

filter_subdist_fb <- function(fname, which.subs, dist.paths, memlimit, 
                              verbose=TRUE, parallel=FALSE, 
                              type="double", overwrite=FALSE, 
                              gower=FALSE, gower.paths=NULL)
{
    vcat(verbose, "Filtering subjects in distance matrices")
    
    if (!is.character(fname))
        stop("input must be a filename")
    if (gower && is.null(gower.paths))
        stop("must supply gower.paths if giving gower")
    
    vcat(verbose, "...reading input")
    input.path <- dirname(fname)
    sdist1 <- attach.big.matrix(fname)
    
    vcat(verbose, "...setup needed variables")
    nr1 <- nrow(sdist1); nobs1 <- sqrt(nr1)
    nc <- ncol(sdist1)
    nobs2 <- length(which.subs); nr2 <- nobs2^2
    if (nobs2 > nobs1)
        stop("length error for which.subs: no filtering can be done")
    matinds <- matrix(1:nr1, nobs1, nobs1)
    matinds <- as.vector(matinds[which.subs,which.subs])
    progress <- ifelse(verbose, "text", "none")
    
    # determine chunks to process
    memlimit <- as.numeric(memlimit)
    vcat(verbose, "...constraining given memory limit of %.2f GB", memlimit)
    
    mem_dmat1 <- n2gb(nr1*1)
    mem_dmat2 <- n2gb(nr2*1)
    f <- function(v) memlimit - v*mem_dmat1 - v*mem_dmat2
    
    if (f(nc) > 0) {
        blocksize <- nc
    } else {
        v <- floor(tryCatch(uniroot(f, c(2,nc))$root, error=function(ex) NA))
        if (is.na(v) || v < 1)
            stop("you don't have enough memory, consider changing --memlimit")
        blocksize <- v
    }
    if (blocksize < 1)
        stop("you do not have enough memory")
    blocks <- niftir.split.indices(1, nc, by=blocksize)
    
    vcat(verbose, "...creating new gower matrices")
    ## setup
    outfile <- file.path(dist.paths$bpath, dist.paths$dfile)
    if (file.exists(outfile) && !overwrite)
        vstop("output '%s' already exists", dist.paths$bpath)
    ## create
    sdist2 <- big.matrix(nr2, nc, type=type, 
                         backingpath=dist.paths$bpath, 
                         backingfile=dist.paths$bfile, 
                         descriptorfile=dist.paths$dfile)
    if (gower) {
        gdist2 <- big.matrix(nr2, nc, type=type, 
                             backingpath=gower.paths$bpath, 
                             backingfile=gower.paths$bfile, 
                             descriptorfile=gower.paths$dfile)
    }
    
    vcat(verbose, "...looping through in %i blocks", blocks$n)
    for (bi in 1:blocks$n) {
        vcat(verbose, "...block %i", bi)
        
        fCol <- blocks$starts[bi]; lCol <- blocks$ends[bi]
        sub_cols <- fCol:lCol; sub_ncol <- length(sub_cols)
        
        # copy
        vcat(verbose, "\t...copying")
        tmp_sdist1 <- deepcopy(x=sdist1, rows=matinds, cols=sub_cols, 
                               type=type, shared=parallel)
        sdist1 <- free.memory(sdist1, input.path)
        
        # save
        vcat(verbose, "\t...saving")
        sub_sdist2 <- sub.big.matrix(sdist2, firstCol=fCol, lastCol=lCol, 
                                     backingpath=dist.paths$bpath)
        deepcopy(x=tmp_sdist1, y=sub_sdist2)
        flush(sub_sdist2); flush(sdist2)
        sdist2 <- free.memory(sdist2, dist.paths$bpath)
        
        # gower?
        if (gower) {
            vcat(verbose, "\t...gowering")
            tmp_gdist1 <- gower.subdist2(tmp_sdist1, verbose=0, parallel=parallel)
            rm(tmp_sdist1); gc(FALSE, TRUE)
            
            # save
            vcat(verbose, "\t...saving gower")
            sub_gdist2 <- sub.big.matrix(gdist2, firstCol=fCol, lastCol=lCol, 
                                         backingpath=gower.paths$bpath)
            deepcopy(x=tmp_gdist1, y=sub_gdist2)
            flush(sub_gdist2); flush(gdist2)
            gdist2 <- free.memory(gdist2, gower.paths$bpath)
            
            rm(tmp_gdist1); gc(FALSE, TRUE)
        } else {
            rm(tmp_sdist1); gc(FALSE, TRUE)
        }
    }
    
    return(sdist2)
}


### NEW CODE USING ARMADILLO

compute_subdist <- function(funclist, subdist, seed_inds, blocksize, ztransform, start=1, verbose=TRUE, testonly=FALSE) {
    nseeds <- length(seed_inds)
    blocks <- niftir.split.indices(start, nseeds, by=blocksize)
    
#    dfun <- function(i, blocks, seed_inds, funclist, subdist, ztransform, verbose, pb) {
    dfun <- function(i, ...) {
        if (verbose) {
            update(pb, i)
            #msg <- sprintf("\nblock %i with voxels %i:%i\n", i, blocks$starts[i], blocks$ends[i])
            #cat(msg)
        }
        dist_inds_CHUNK <- blocks$starts[i]:blocks$ends[i]
        scor_inds_CHUNK <- seed_inds[dist_inds_CHUNK]
        cormaps_list <- vbca_batch(funclist, scor_inds_CHUNK, ztransform=ztransform, shared=FALSE)        
        tmp <- compute_subdist_worker(cormaps_list, scor_inds_CHUNK, subdist, dist_inds_CHUNK)
        rm(dist_inds_CHUNK, scor_inds_CHUNK, cormaps_list, tmp)
        gc(FALSE, TRUE)
        return(NULL)
    }
    
    # Test
    i <- 1
    if (verbose) {
        cat("...running a test (", blocks$starts[i],  ")\n")
        pb <- progressbar(i)
    } else {
        pb <- NULL
    }
    dfun(i)
    check_dmat(matrix(subdist[,blocks$starts[i]], sqrt(nrow(subdist))))
    check_dmat(matrix(subdist[,blocks$ends[i]], sqrt(nrow(subdist))))
    if (verbose)
        end(pb)
    if (testonly) {
        cat("...test only...\n")
        return(NULL)
    }
    
    # Subdist Calculation
    if (verbose) {
        cat("...now the real deal\n")
        pb <- progressbar(blocks$n)
    } else {
        pb <- NULL
    }
    
    if (getDoParRegistered() && getDoParWorkers() > 1) {
        lo <- min(getDoParWorkers()*3, blocks$n-1)
        superblocks <- niftir.split.indices(2, blocks$n, length.out=lo)
        foreach(si=1:superblocks$n, .packages=c("connectir"), .inorder=TRUE) %dopar% 
            for(i in superblocks$starts[si]:superblocks$ends[si]) dfun(i)
    }
    else {
        for (i in 2:blocks$n)
            dfun(i)
    }
    
    if (verbose)
        end(pb)
}

compute_subdist_worker <- function(sub.cormaps, cor_inds, outmat, dist_inds, type="double", ...) {
    nsubs <- length(sub.cormaps)
    nvoxs <- ncol(sub.cormaps[[1]])
    nseeds <- nrow(sub.cormaps[[1]])
    if (nseeds != length(cor_inds) || nseeds != length(dist_inds))
        stop("length of inds doesn't match nrow of first sub.cormaps element")
    
    #if (is.null(outmat))
    #    outmat <- big.matrix(nsubs^2, nseeds, type=type, ...)
    #else if (ncol(outmat) != nseeds || nrow(outmat) != nsubs^2)
    #    stop("dimensions of outmat do not match nsubs and nseeds values")
    
    subsMap <- big.matrix(nvoxs-1, nsubs, type=type, shared=FALSE, ...)
    ALPHA <- 1/(nvoxs-2)
    voxs <- 1:nvoxs
    for (i in 1:nseeds) {
        .Call("CombineSubMapsMain", sub.cormaps, subsMap@address, as.double(i), as.double(voxs[-cor_inds[i]]), as.double(nvoxs-1), as.double(nsubs))
        col <- sub.big.matrix(outmat, firstCol=dist_inds[i], lastCol=dist_inds[i])
        dgemm(C=col, A=subsMap, B=subsMap, TRANSA='t', ALPHA=ALPHA, LDC=as.double(nsubs))
        .Call("BigSubtractScalarMain", col@address, as.double(1), TRUE);
    }
    
    rm(subsMap)
    gc(F)
    
    return(outmat)
}
