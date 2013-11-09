# Higher-level wrappers around certain R functions

roi_mean_wrapper <- function(func_file, roi_file, mask_file=NULL, 
                              out_file=NULL, outtype="nifti", 
                              to_return=FALSE, overwrite=FALSE, 
                              verbose=TRUE, shared=FALSE)
{
    vcat(verbose, "Averaging mean signal in '%s' with '%s' ROIs", 
         basename(func_file), basename(roi_file))
    
    if (!is.character(func_file) || !file.exists(func_file))
        stop("Could not find functional file ", func_file)
    if (!is.character(roi_file) || !file.exists(roi_file))
        stop("Could not find roi file ", roi_file)
    if (!is.null(out_file) && file.exists(out_file) && !overwrite) {
        vcat(verbose, "File '%s' already exists, not re-running", out_file)
        return(NULL)
    }
    if (is.null(out_file) && !to_return)
        stop("You haven't specified an output file or asked to return the data")
    if (!(outtype %in% c("nifti", "text")))
        stop("output type can only be nifti or text")
    
    vcat(verbose, "...reading data")
    func <- read.big.nifti4d(func_file, shared=shared)
    rois <- read.mask(roi_file, NULL)
    hdr <- read.nifti.header(func_file)
    if (!is.null(mask_file)) {
        mask <- read.mask(mask_file) & rois!=0
    } else {
        mask <- rois!=0
    }
    
    vcat(verbose, "...masking")
    func <- do.mask(func, mask, shared=shared)
    rois <- rois[mask]
    
    vcat(verbose, "...averaging")
    new_func <- roi_mean(func, rois)
    
    if (!is.null(out_file)) {
        vcat(verbose, "...writing file '%s'", out_file)
        if (outtype == "nifti") {
            hdr$dim <- dim(new_func)
            hdr$pixdim <- c(hdr$pixdim[4], 1)
            write.nifti(new_func, hdr, odt="float", 
                        outfile=out_file, overwrite=overwrite)
        } else if (outtype == "text") {
            write.table(new_func, file=out_file, quote=FALSE, 
                        row.names=F, col.names=unique(rois))
        }
    }
    
    if (to_return) {
        vcat(verbose, "...returning data")
        return(new_func)
    }
}

wrap_gcor <- function(func_file, mask_file, out_file=NULL, 
                      blocksize=0, memlimit=4, 
                      to_return=FALSE, overwrite=FALSE, 
                      verbose=TRUE, parallel=FALSE, shared=parallel, 
                        ...) 
{
    vcat(verbose, "Running global correlation for '%s'", func_file)
    
    if (!is.character(func_file) || !file.exists(func_file))
        stop("Could not find functional file ", func_file)
    if (!is.character(mask_file) || !file.exists(mask_file))
        stop("Could not find mask file ", mask_file)
    if (!is.null(out_file) && file.exists(out_file) && !overwrite) {
        vcat(verbose, "File '%s' already exists, not re-running", out_file)
        return(NULL)
    }
    
    vcat(verbose, "...reading mask")
    mask <- read.mask(mask_file)
    hdr <- read.nifti.header(mask_file)
    
    vcat(verbose, "...calculating memory demands")
    tmp <- read.nifti.header(func_file)
    if (length(tmp$dim) != 4)
        vstop("functional file '%s' must be 4D", func_file)
    blocksize <- get_gcor_limit(blocksize, memlimit, sum(mask), tmp$dim[[4]], verbose)
    
    vcat(verbose, "...reading functional")
    bmat <- load_and_mask_func_data(func_file, mask, type="double", shared=shared)
    
    vcat(verbose, "...correlating")
    vec <- gcor(bmat, blocksize, verbosity=verbose*1, parallel=parallel, shared=shared, 
                ...)
    
    if (!is.null(out_file)) {
        vcat(verbose, "...writing file '%s'", out_file)
        write.nifti(vec, hdr, mask, outfile=out_file, overwrite=overwrite)
    }
    
    if (to_return) {
        vcat(verbose, "...returning data")
        return(vec)
    }
}


wrap_reho <- function(func_file, mask_file, out_file=NULL, 
                      to_return=FALSE, overwrite=FALSE, 
                      verbose=TRUE, parallel=FALSE, shared=parallel, 
                      ...) 
{
    vcat(verbose, "Running regional homogeneity for '%s'", func_file)
    
    if (!is.character(func_file) || !file.exists(func_file))
        stop("Could not find functional file ", func_file)
    if (!is.character(mask_file) || !file.exists(mask_file))
        stop("Could not find mask file ", mask_file)
    if (!is.null(out_file) && file.exists(out_file) && !overwrite) {
        vcat(verbose, "File '%s' already exists, not re-running", out_file)
        return(NULL)
    }
    
    vcat(verbose, "...reading mask")
    mask <- read.mask(mask_file)
    hdr <- read.nifti.header(mask_file)
    
    vcat(verbose, "...reading functional")
    bmat <- load_and_mask_func_data3(func_file, mask, type="double", shared=shared)
    
    vcat(verbose, "...rehoing")
    res <- reho(bmat, verbose=verbose, parallel=parallel, ...)
    
    if (!is.null(out_file)) {
        vcat(verbose, "...writing file '%s'", out_file)
        write.nifti(res$reho, res$hdr, res$mask, outfile=out_file, overwrite=overwrite)
    }
    
    if (to_return) {
        vcat(verbose, "...returning data")
        return(res)
    }
}

.wrap_functionals <- function(func_files1, mask_file1, 
                              func_files2=NULL, mask_file2=NULL, 
                              verbose=TRUE, parallel=FALSE, shared=parallel, 
                              extra_checks=FALSE, 
                              memfunc=NULL, memlimit=1, ...)
{
    vcat(verbose, "Reading inputs")
    
    if (!is.null(mask_file2) && is.null(func_files2))
        func_files2 <- func_files1
    .check_mask_paths(mask_file1, mask_file2)
    .check_func_paths(func_files1, func_files2)
    
    if (!is.null(mask_file1))
        hdr <- read.nifti.header(mask_file1)
    else
        hdr <- read.nifti.header(func_files1[[1]])
    
    vcat(verbose, "...reading mask")
    mask1 <- .get_mask(func_files1, mask_file1)
    nvoxs1 <- sum(mask1)
    if (!is.null(func_files2)) {
        mask2 <- .get_mask(func_files2, mask_file2)
        nvoxs2 <- sum(mask2)
    } else {
        mask2 <- NULL
        nvoxs2 <- NULL
    }
        
    vcat(verbose, "...getting # of time-points for functional data")
    ntpts <- .get_nifti_nvols(func_files1, func_files2, verbose)
    
    if (!is.null(memlimit) || !is.null(memfunc)) {
        vcat(verbose, "...calculating memory demands")
        blocksize <- memfunc(0, memlimit, nvoxs1, ntpts, nvoxs2, verbose, ...)
    } else {
        blocksize <- NULL
    }
    
    vcat(verbose, "...reading %i functionals", length(func_files1))
    ret <- .read_funcs(func_files1, mask1, func_files2, mask2, 
                       verbose, parallel, extra_checks)
    
    ret$blocksize <- blocksize
    ret$hdr <- hdr
    ret$mask <- mask1
    
    return(ret)
}

wrap_functionals <- function(...) .wrap_functionals(...)

wrap_kendall <- function(func_files1, mask_file1, 
                            func_files2=NULL, mask_file2=NULL, 
                            design_mat=NULL, 
                            to_return=FALSE, out_file=NULL, overwrite=FALSE, 
                            verbose=TRUE, parallel=FALSE, shared=parallel, 
                            memlimit=4, extra_checks=FALSE, ...)
{
    vcat(verbose, "Running kendall's W")
    progress <- ifelse(verbose, "text", "none")
    
    if (is.null(out_file) && !to_return)
        stop("Must specify at least one of 'out_file' or 'to_return'")
    if (!is.null(out_file) && file.exists(out_file) && !overwrite) {
        vcat(verbose, "File '%s' already exists, not re-running", out_file)
        return(NULL)
    }
    
    ret <- .wrap_functionals(func_files1, mask_file1, 
                             func_files2, mask_file2, 
                             verbose, parallel, shared, 
                             extra_checks, 
                             .get_kendall_limit, memlimit)
    
    if (!is.null(design_mat)) {
        vcat(verbose, "reading in design matrix")
        tmp_fname <- design_mat
        tmp <- as.matrix(read.table(design_mat, header=TRUE))
        design_mat <- big.matrix(nrow(tmp), ncol(tmp), type="double", shared=TRUE)
        design_mat[,] <- tmp[,]; rm(tmp)
        k <- qlm_rank(design_mat)
        if (k < ncol(design_mat))
            vstop("design matrix (--regress %s) is rank deficient", tmp_fname)
        rm(tmp_fname)
    }
    
    vcat(verbose, "\nKendalling")
    verbosity <- ifelse(verbose, 2, 0)
    start.time <- Sys.time()
    vec <- kendall3(ret$funcs1, ret$blocksize, ret$funcs2, design_mat, 
                   verbose=verbose, parallel=parallel, ...)
    end.time <- Sys.time()
    vcat(verbose, "Kendalling is done! It took: %.2f minutes\n", 
         as.numeric(end.time-start.time, units="mins"))
    
    if (!is.null(out_file)) {
        vcat(verbose, "...writing file '%s'", out_file)
        write.nifti(vec, ret$hdr, ret$mask, outfile=out_file, 
                    overwrite=overwrite, odt="float")
    }
    
    if (to_return) {
        vcat(verbose, "...returning data")
        return(vec)
    }
}

wrap_glm <- function(func_files1, mask_file1, ev_file, contrast_file, 
                     outdir, summarize=FALSE, overwrite=FALSE, 
                     func_files2=NULL, mask_file2=NULL, 
                     memlimit=4, extra_checks=FALSE, 
                     verbose=TRUE, parallel=FALSE, shared=parallel, 
                     ztransform=FALSE) 
{
    vcat(verbose, "Running GLM")
    progress <- ifelse(verbose, "text", "none")
    
    vcat(verbose, "...reading and setting up EV and contrast files")
    if (!is.character(ev_file) || !file.exists(ev_file)) 
        stop("Could not find EV file ", ev_file)
    if (!is.character(contrast_file) || !file.exists(contrast_file)) 
        stop("Could not find contrast file ", contrast_file)
    evs <- as.matrix(read.table(ev_file, header=T))
    cons <- as.matrix(read.table(contrast_file, header=T))
    contrast_names <- rownames(cons)
    if (is.numeric(contrast_names))
        stop("must have row names for contrast matrix")
    if (ncol(evs) != ncol(cons)) {
        vstop(paste("# of columns in EV file '%s' not the same as # of columns",  
                    "in contrast file '%s'"), ev_file, contrast_file)
    }
    tmp <- big.matrix(nrow(evs), ncol(evs), type="double", shared=shared)
    tmp[,] <- evs; evs <- tmp; rm(tmp)
    k <- qlm_rank(evs)
    if (k < ncol(evs))
        vstop("EV file '%s' is rank deficient", ev_file)
    nevs <- ncol(evs)
    ncons <- nrow(cons)
    
    # Read in functionals
    ret <- wrap_functionals(func_files1, mask_file1, 
                             func_files2, mask_file2, 
                             verbose, parallel, shared, 
                             extra_checks, 
                             .get_glm_limit, memlimit, nevs, ncons)
    
    vcat(verbose, "...setting up output")
    outdir <- abspath(outdir)
    if (!is.null(outdir) && file.exists(outdir)) {
        if (!overwrite) {
            vcat(verbose, "Directory '%s' already exists, not re-running", outdir)
            return(NULL)
        } else {
            vcat(verbose, "Trying to remove directory '%s'", outdir)
            vsystem("rm %s/infuncs/*", outdir)
            vsystem("rmdir %s/infuncs", outdir)
            vsystem("rm %s/summary/*", outdir)
            vsystem("rmdir %s/summary", outdir)
            vsystem("rm %s/*", outdir)
            vsystem("rmdir %s", outdir)
        }
    } else {
        dir.create(outdir)
    }
    ## copy evs
    file.copy(ev_file, file.path(outdir, "model_evs.txt"))
    ## copy contrasts
    file.copy(contrast_file, file.path(outdir, "model_contrasts.txt"))
    ## copy mask
    if (!is.null(mask_file1))
        file.copy(mask_file1, file.path(outdir, "mask.nii.gz"))
    #else {
    #    hdr <- read.nifti.header(func_files1[1])
    #    hdr$dim <- hdr$dim[-length(hdr$dim)]; hdr$dim <- hdr$pixdim[-length(hdr$pixdim)]
    #    if (length(hdr$dim) == 1) {
    #        hdr$dim <- c(hdr$dim, 1, 1)
    #        hdr$pixdim <- c(hdr$pixdim, 1, 1)
    #    }
    #    mask <- as.nifti(array(1, hdr$dim), hdr)
    #    write.nifti(mask, odt="int", outfile=file.path(outdir, "mask.nii.gz"))
    #}
    if (!is.null(mask_file2))
        file.copy(mask_file2, file.path(outdir, "mask2.nii.gz"))
    ## soft-link functionals
    dir.create(file.path(outdir, "infuncs"))
    for (i in 1:length(func_files1)) {
        file.symlink(func_files1[i], file.path(outdir, "infuncs", 
                                           sprintf("func%04i.nii.gz", i)))
    }
    ## summary output (only given if summarize=TRUE)
    dir.create(file.path(outdir, "summary"))
    
    # Do It!
    vcat(verbose, "...GLMing")
    start.time <- Sys.time()
    outmats <- vox_glm3(ret$funcs1, evs, cons, ret$blocksize, ret$funcs2, bp=outdir, 
                         verbose=verbose, parallel=parallel, shared=shared, 
                         ztransform=ztransform)
    end.time <- Sys.time()
    vcat(verbose, "GLMing is done! It took: %.2f minutes\n", 
         as.numeric(end.time-start.time, units="mins"))
    
    # Summarize
    df <- nrow(evs) - ncol(evs)
    if (df < 1)
        stop("non-positive degrees of freedom")
    if (summarize) {
        vcat(verbose,  "...summarizing results")
        start.time <- Sys.time()
        sdir <- file.path(outdir, "summary")
        glm_summarize(tmats=outmats, tdir=outdir, df=df,
                      outdir=sdir, outhdr=ret$hdr, outmask=ret$mask, 
                      memlimit=memlimit, 
                      verbose=verbose, parallel=parallel, shared=shared)
        end.time <- Sys.time()
        vcat(verbose, "Summarizing is done! It took: %.2f minutes\n", 
             as.numeric(end.time-start.time, units="mins"))
    }
    
    invisible(outmats)
}

wrap_glmnet_subdist_cross <- function(sdist_file, mask_file, label_file, 
                                    out_prefix=NULL, to_return=FALSE, overwrite=FALSE, 
                                    family=NULL, standardize=TRUE, cross=10, 
                                    memlimit=1, verbose=TRUE, parallel=FALSE)
{
    if (!file.exists(sdist_file))
        vstop("Subject Distances '%s' does not exist", sdist_file)
    if (!file.exists(mask_file))
        vstop("Mask file '%s' does not exist", mask_file)
    if (!file.exists(label_file))
        vstop("Label file '%s' does not exist", label_file)
    if (is.null(out_prefix) && !to_return)
        vstop("Must specify either out_file or to_return")
    
    if (!is.null(out_prefix)) {
        out_files <- c("error_min", "error_mean", "nzero_min", "nzero_mean")
        out_files <- sprintf("%s_%s.nii.gz", out_prefix, out_files)
        for (out_file in out_files) {
            if (!is.null(out_file) && file.exists(out_file) && !overwrite)
                vstop("Output file '%s' already exists (must specify overwrite)", out_file)
        }
    }
    
    bpath <- dirname(sdist_file)
    Xs <- attach.big.matrix(sdist_file)
    y <- read.table(label_file)[,1]
    
    mask <- read.mask(mask_file)
    hdr <- read.nifti.header(mask_file)
    
    res <- glmnet_subdist_cross(Xs, y, family=family, standardize=standardize, 
                                cross=cross, memlimit=memlimit, verbose=verbose, 
                                parallel=parallel, bpath=bpath)
    res <- sapply(res, function(x) x)
    save(res, file="tmp.rda")
    
    if (!is.null(out_prefix)) {
        vcat(verbose, "...saving")
        for (i in 1:length(out_files))
            write.nifti(res[,i], hdr, mask, odt="float", outfile=out_files[i])
    }
    
    if (to_return)
        return(res)
}

wrap_svm_subdist_cross <- function(sdist_file, mask_file, label_file, 
                                    out_file=NULL, to_return=FALSE, overwrite=FALSE, 
                                    kernel="linear", type=NULL, cross=10, 
                                    memlimit=1, verbose=TRUE, parallel=FALSE)
{
    if (!file.exists(sdist_file))
        vstop("Subject Distances '%s' does not exist", sdist_file)
    if (!file.exists(mask_file))
        vstop("Mask file '%s' does not exist", mask_file)
    if (!file.exists(label_file))
        vstop("Label file '%s' does not exist", label_file)
    if (is.null(out_file) && !to_return)
        vstop("Must specify either out_file or to_return")
    if (!is.null(out_file) && file.exists(out_file) && !overwrite)
        vstop("Output file '%s' already exists (must specify overwrite)", out_file)
    
    bpath <- dirname(sdist_file)
    Xs <- attach.big.matrix(sdist_file)
    y <- read.table(label_file)[,1]
    
    mask <- read.mask(mask_file)
    hdr <- read.nifti.header(mask_file)
    
    if (is.null(type)) {
        res <- svm_subdist_cross(Xs, y, kernel=kernel, cross=cross, memlimit=memlimit, verbose=verbose, parallel=parallel, bpath=bpath)
    } else {
        res <- svm_subdist_cross(Xs, y, kernel=kernel, cross=cross, memlimit=memlimit, verbose=verbose, parallel=parallel, bpath=bpath, type=type)
    }
    
    if (!is.null(out_file)) {
        vcat(verbose, "...saving")
        write.nifti(res, hdr, mask, odt="float", outfile=out_file)
    }
    
    if (to_return)
        return(res)
}

wrap_kmeans_subdist_cross <- function(sdist_file, mask_file, label_file, 
                                      out_file=NULL, to_return=FALSE, overwrite=FALSE, 
                                      iter.max=200, nstart=20, algorithm="Hartigan-Wong", 
                                      memlimit=1, verbose=TRUE, parallel=FALSE)
{
    if (!file.exists(sdist_file))
        vstop("Subject Distances '%s' does not exist", sdist_file)
    if (!file.exists(mask_file))
        vstop("Mask file '%s' does not exist", mask_file)
    if (!file.exists(label_file))
        vstop("Label file '%s' does not exist", label_file)
    if (is.null(out_file) && !to_return)
        vstop("Must specify either out_file or to_return")
    if (!is.null(out_file) && file.exists(out_file) && !overwrite)
        vstop("Output file '%s' already exists (must specify overwrite)", out_file)
    
    bpath <- dirname(sdist_file)
    Xs <- attach.big.matrix(sdist_file)
    y <- as.numeric(read.table(label_file)[,1])
    
    mask <- read.mask(mask_file)
    hdr <- read.nifti.header(mask_file)
    
    res <- kmeans_subdist_cross(Xs, y, iter.max=iter.max, nstart=nstart, 
                                 algorithm=algorithm, memlimit=memlimit, 
                                 verbose=verbose, parallel=parallel, bpath=bpath)
    
    if (!is.null(out_file)) {
        vcat(verbose, "...saving")
        write.nifti(res, hdr, mask, odt="float", outfile=out_file)
    }
    
    if (to_return)
        return(res)
}
