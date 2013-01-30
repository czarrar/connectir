# Wrapper read functions
read.big.tab <- function(file, ...) read.big.matrix(file, sep='\t', ...)

read.big.space <- function(file, ...) read.big.matrix(file, sep=' ', ...)

read.big.csv <- function(file, ...) read.big.matrix(file, sep=',', ...)

read.big.matlab <- function(file, ...) {
    library(R.matlab)
    m <- readMat(file)
    if (length(m) > 1)
        stop(sprintf("matlab input '%s' can't have more than 1 variable, but %i found", 
                        fname, length(m)))
    as.big.matrix(m[[1]], ...)
}

# assume that rows are time-points & columns are nodes/voxels
read.big.nifti2d <- function(file, ...) {
    img <- read.nifti.image(file)
    if (length(dim(img)) == 1)
        dim(img) <- c(length(img), 1)
    bigmat <- as.big.matrix(img, ...)
    return(bigmat)
}

# Read input data as big matrix
gen_big_reader <- function(intype, ...) {
    choices <- c("nifti2d", "nifti4d", "nifti", "tab", "space", "csv", "matlab")
    if (!(intype %in% choices))
        stop("unrecognized input type: ", intype)
    
    args <- list(...)
    fun <- function(x, ...) {
        # Checks
        if (is.big.matrix(x))
            return(x)
        else if (!is.character(x))
            stop("input must be a character or big.matrix and not ", class(x))
        
        # Read
        args$file <- x
        mat <- do.call(sprintf("read.big.%s", intype), args)
        
        gc(FALSE, TRUE)
        return(mat)
    }
    
    return(fun)
}

# Detect the type of files
detect_ftypes <- function(fnames, force.type=NULL, verbose=TRUE) {
    df <- data.frame(
        extensions = c(".nii.gz", ".nii", ".hdr", ".img", 
                       ".mat", ".txt", ".1D", ".tab", ".csv"), 
        formats = c("nifti", "nifti", "nifti", "nifti", 
                    "matlab", "space", "space", "tab", "csv")
    )
    
    if (!is.null(force.type) && !(force.type %in% df$formats))
        stop("unrecognized file type ", force.type)
    
    formats <- sapply(fnames, function(fname) {
        find <- sapply(df$extensions, function(ext) grepl(ext, fname))
        if (any(find))
            return(as.character(df$formats[find]))
        else if (is.null(force.type))
            stop("unrecognized extension in ", fname)
    })
    
    if (any(formats != formats[1]))
        warning("not all extensions are ", formats[1], " like in ", fnames[1], immediate.=TRUE)
    if (!is.null(force.type) && any(formats != force.type))
        vcat(verbose, "not all extensions are the same as the forced one ", force.type)
    
    ftype <- ifelse(is.null(force.type), formats[1], force.type)
    
    return(ftype)
}

# Only get time-points for set of functional files
get_funclist_tpts <- function(inlist) {
    sapply(1:length(inlist$files), function(i) {
        if (inlist$ftype == "nifti") {
            hdr <- read.nifti.header(inlist$files[i])
            return(hdr$dim[4])
        } else if (inlist$ftype %in% c("space", "tab", "csv")) {
            cmd <- sprintf("wc -l %s | awk '{print $1}'", inlist$files[i])
            return(system(cmd, intern=T))
        } else {
            vstop("Cannot read # of time-points for type '%s'", inlist$ftype)
        }
            
    })
}

mask.read <- function(file, ftype)
{
    if (ftype == "nifti") {
        mask <- read.nifti.image(file)
        mask <- as.vector(mask)
    } else if (ftype == "space") {
        mask <- read.table(file)
        mask <- as.vector(as.matrix(mask))
    } else {
        vstop("unrecognized file type '%s' for mask '%s'", ftype, file)
    }
    mask > 0
}

mask.auto <- function(file, read_fun, na.rm=FALSE, ...)
{
    x <- read_fun(file, ...)
    mask <- colvar(x, na.rm=na.rm)
    mask != 0
}

# If automask is false then files are assumed to be 4D or 2D
# If automask is true then files are assumed to be 3D or 1D
overlap_masks <- function(files, automask=FALSE, verbose=TRUE)
{
    vcat(verbose, "Calculating overlap of masks")
    
    progress <- ifelse(verbose, "text", "none")
    
    vcat(verbose, "...detecting file type")
    ftype <- detect_ftypes(files, verbose=verbose)
    
    if (automask) {
        vcat(verbose, "...reading in functionals and automasking")
        reader <- gen_big_reader(ftype, type="double", shared=FALSE)
        # ghetto fix for not setting the type properly...todo
        fun <- function(...) suppressWarnings(mask.auto(...))
        masks <- laply(files, fun, reader, .progress=progress)
    } else {
        vcat(verbose, "...reading in masks")
        masks <- laply(files, mask.read, ftype, .progress=progress)
    }
    
    vcat(verbose, "...checking for NaNs")
    nas <- is.na(masks)
    if (any(nas))
        vstop("%i NaNs found", sum(nas))
    
    vcat(verbose, "...calculating overlap")
    overlap <- colMeans(masks)
    
    vcat(verbose, "...creating overlap mask")
    mask <- overlap == 1
    percent_nonzero <- round(mean(mask)*100, 2)
    vcat(verbose, "mask has %.2f%% of non-zero voxels", percent_nonzero)
    
    return(list(mask=mask, overlap=overlap, sub.masks=masks))
}




exclude_sub_masks <- function(inmasks, exclude.thresh=0.05, filter.mask=NULL, 
                              verbose=FALSE, parallel=FALSE) 
{
    tmp <- read.mask(inmasks[[1]])
    nsubs <- length(inmasks)
    nnodes <- length(tmp)
    progress <- ifelse(verbose, "text", "none")
    
    all.masks <- big.matrix(nsubs, nnodes)
    
    vcat(verbose, "reading in data")
    tmp <- llply(1:nsubs, function(i) {
        all.masks[i,] <- read.mask(inmasks[[i]])
        return(NULL)
    }, .progress=progress, .parallel=parallel)
    
    vcat(verbose, "computing overlap")
    overlap <- colmean(all.masks)
    
    vcat(verbose, "filtering")
    if (!is.null(filter.mask)) {
        filter <- read.mask(filter.mask) & overlap > 0
    } else {
        filter <- overlap > 0
    }
    nnodes <- sum(filter)
    overlap <- overlap[filter]
    all.masks <- deepcopy(all.masks, cols=which(filter))
    
    vcat(verbose, "getting bad subjects to exclude")
    sub.ns <- big.matrix(nsubs, 1, init=0)
    tmp <- llply(1:nnodes, function(i) {
        w <- which(all.masks[,i]!=1)
        sub.ns[w,1] <- sub.ns[w,1] + 1
        return(NULL)
    }, .progress=progress, .parallel=parallel)
    sub.ns <- sub.ns[,1]/nnodes
    exclude.subs <- which(sub.ns > exclude.thresh)
    include.subs <- which(sub.ns < exclude.thresh)
    
    vcat(verbose, "re-calculating overlap with only included subjects")
    use.all.masks <- deepcopy(all.masks, rows=include.subs)
    overlap2 <- colmean(use.all.masks)
    
    vcat(verbose, "getting mask")
    mask <- overlap2 == 1
    vcat(verbose, "mask has %i good nodes and %i bad ones", sum(mask), sum(!mask))
    
    list(mask=mask, overlap=overlap2, all=all.masks, 
         bad.ns=sub.ns/nnodes, include=include.subs, exclude=exclude.subs)
}


# Loading Functional Data

## Prepares the way forward
## 1. What are the type of files?
## 2. What is the function to read these files?

## A bunch of mask related stuff
## 3. Do we need to automask the data?
## 4. Do we get the overlap across subject masks?
## 5. Is there a group mask to consider as well?

## 6. Actually load and mask the data

## 6. really should be for any vector => voxelwise transform
## 6. For ROI data, is there a voxelwise ROI map? This should be unrelated to above.


## 4. Read in the functional data given a brain mask


## Prepares the way forward
## 1. What are the type of files?
## 2. What is the function to read these files?
load_funcs.prepare <- function(files, verbose=TRUE) 
{
    vcat(verbose, "Preparing functional info")
    
    # What is the file type?
    ftype <- detect_ftypes(files, verbose=verbose)
    
    # What is the function to read these functional files?
    reader <- gen_big_reader(ftype, type="double")
    
    ret <- list(files=files, ftype=ftype, reader=reader)
    class(ret) <- "inlist_prepare"
    
    return(ret)
}

## A bunch of mask related stuff
## 3. Do we need to automask the data?
## 4. Do we get the overlap across subject masks?
## 5. Is there a group mask to consider as well?
load_funcs.mask <- function(inlist, automask=FALSE, subject.masks=NULL,  
                            group.mask=NULL, verbose=TRUE, detailed=FALSE)
{
    vcat(verbose, "Computing brain mask")
    
    if (class(inlist) != "inlist_prepare")
        stop("Input inlist must be of class inlist_prepare")
    
    mat <- inlist$reader(inlist$files[1])
    default <- rep(T, ncol(mat))
    mask <- auto <- subjects <- group <- default
    rm(mat)
    
    if (automask) {
        auto <- overlap_masks(inlist$files, automask, verbose)
        mask <- mask & auto$mask
    }
    
    if (!is.null(subject.masks)) {
        subjects <- overlap_masks(subject.masks, verbose=verbose)
        mask <- mask & subjects$mask
    }
    
    if (!is.null(group.mask)) {
        ftype <- detect_ftypes(group.mask, verbose=FALSE)
        group <- mask.read(group.mask, ftype)
        percent_nonzero <- round(mean(group)*100, 2)
        mask <- mask & group
        vcat(verbose, "group mask has %.2f%% non-zero voxels", percent_nonzero)
    }
    
    percent_nonzero <- round(mean(mask) * 100, 2)
    vcat(verbose, "final mask has %.2f%% non-zero voxels", percent_nonzero)
    
    if (sum(mask) == 0)
        stop("no non-zero voxels found in final mask")
    
    inlist$mask <- mask
    if (detailed)
        inlist$mask_details <- list(auto=auto, subjects=subjects, group=group)
    class(inlist) <- "inlist"
    
    return(inlist)
}

load_funcs.read <- function(inlist, verbose=TRUE, ...)
{
    vcat(verbose, "Reading in functional data")
    
    if (class(inlist) != "inlist")
        stop("input argument must be of class inlist")
    progress <- ifelse(verbose, "text", "none")
    
    inds <- which(inlist$mask)
    inlist$funcs <- llply(inlist$files, function(f) {
        if (inlist$ftype == "nifti") {
            func <- inlist$reader(f, voxs=inds, ...)
        } else {
            func <- inlist$reader(f, ...)
            func <- deepcopy(func, cols=inds)
        }
    }, .progress=progress)
    
    return(inlist)
}

load_funcs.scale <- function(inlist, verbose=TRUE, to.copy=FALSE, 
                             parallel=FALSE, ...)
{
    vcat(verbose, "Scaling functional data")
    
    if (class(inlist) != "inlist")
        stop("input argument must be of class inlist")
    
    progress <- ifelse(verbose, "text", "none")
    
    if (to.copy)
        inlist$funcs_nonscaled <- inlist$funcs
    
    inlist$funcs <- llply(inlist$funcs, function(func) {
        scale(func, to.copy=to.copy, ...)
    }, .progress=progress, .parallel=parallel)
    
    invisible(gc(FALSE, TRUE))
    
    return(inlist)
}

load_funcs.read_and_scale <- function(inlist, verbose=TRUE, to.copy=FALSE, 
                                      parallel=FALSE, ...)
{
    inlist <- load_funcs.read(inlist, verbose, ...)
    inlist <- load_funcs.scale(inlist, verbose, to.copy, parallel, ...)
    return(inlist)
}

load_funcs <- function(files, verbose=TRUE, 
                       automask=FALSE, subject.masks=NULL,  
                       group.mask=NULL, detailed=FALSE, 
                       to.copy=FALSE, parallel=FALSE, 
                       ...)
{
    inlist <- load_funcs.prepare(files, verbose)
    inlist <- load_funcs.mask(inlist, automask, subject.masks,  
                              group.mask, verbose, detailed)
    inlist <- load_funcs.read(inlist, verbose, ...)
    inlist <- load_funcs.scale(inlist, verbose, to.copy, parallel, ...)
    return(inlist)             
}


# Load data
load_and_mask_func_data2 <- function(xs, read_fun, mask=NULL, verbose=FALSE, 
                                     scale=TRUE, ...)
{
    if (!is.list(xs) && !is.vector(xs))
        stop("input 'xs' must be a vector or list")
    if (!is.function(read_fun))
        stop("input 'read_fun' must be a function")
    progress <- ifelse(verbose, "text", "none")
    
    if (!is.null(mask) && is.logical(mask))
        mask <- which(mask)
    
    vcat(verbose, "reading data")
    dat.list <- llply(xs, function(x) {
        x <- read_fun(x, ...)
        if (!is.null(mask)) {
            z <- deepcopy(x, cols=mask, ...)
            rm(x); gc(FALSE, TRUE)
            x <- z
            rm(z); gc(FALSE, TRUE)
        }
        y <- scale(x, scale=scale, to.copy=TRUE, ...)
        rm(x); gc(FALSE, TRUE)
        return(y)
    }, .progress=progress, .parallel=FALSE)
    
    gc(FALSE, TRUE)
    
    return(dat.list)
}

check_func_data <- function(xs, dat.list, extra=FALSE, verbose=FALSE, parallel=FALSE) 
{
    vcat(verbose, "checking data")
    progress <- ifelse(verbose, "text", "none")
    nc <- ncol(dat.list[[1]])
    n <- length(dat.list)
    rets <- llply(1:n, function(i) {
        dat <- dat.list[[i]]
        fname <- xs[[i]]
        
        # Everything must have same # of voxels
        if (nc != ncol(dat)) {
            vcat(verbose, 
                 "%s must have the same # of nodes (%i vs %i) as other datasets", 
                  fname, nc, ncol(dat))
            return(1)
        }
        
        if (extra) {
            # There can't be any NaNs
            col.nas <- colna(dat.list[[i]])>0
            if (any(col.nas)) {
                w <- paste(which(col.nas), collapse=", ")
                vcat(verbose, "%s has NaNs in nodes %s", fname, w)
                return(2)
            }
            
            # Extra check (ensure standard deviation greater than 0)
            col.sds <- colsd(dat.list[[1]])==0
            if (any(col.sds)) {
                w <- paste(which(col.sds), collapse=", ")
                vcat(verbose, "%s has 0 sd in nodes %s", fname, w)
                return(3)
            }
        } else {
            # Check for NAs at least for the second row
            nas <- is.na(dat.list[[i]][2,])
            if (any(nas)) {
                w <- paste(which(nas), collapse=", ")
                vcat(verbose, "%s has NaNs in nodes %s", fname, w)
                return(2)
            }
        }
        
        return(0)
    }, .progress=progress, .parallel=parallel)
    
    unlist(rets)
}

###
# HELPER FUNCTIONS
###

.check_mask_paths <- function(mask_file1=NULL, mask_file2=NULL) {
    if (!is.null(mask_file1) && (!is.character(mask_file1) || !file.exists(mask_file1)))
        stop("Could not find mask file ", mask_file1)
    if (!is.null(mask_file2) && (!is.character(mask_file2) || !file.exists(mask_file2)))
        stop("Could not find mask file ", mask_file2)
}

.check_func_paths <- function(func_files1, func_files2=NULL) {
    if (!is.character(func_files1) && length(func_files1) < 2)
        stop("func_files must be a vector of at least 2 filenames")
    if (!is.null(func_files2) && length(func_files1) != length(func_files2))
        stop("1st set of functional files is a different length than 2nd set")
    for (func_file in func_files1) {
        if (!is.character(func_file) || !file.exists(func_file))
            stop("Could not find functional file ", func_file)
    }
    if (!is.null(func_files2)) {
        for (func_file in func_files2) {
            if (!is.character(func_file) || !file.exists(func_file))
                stop("Could not find functional file ", func_file)
        }
    }
}

.get_mask <- function(infiles, mask=NULL) {
    # Check nvoxs in functional
    hdr <- read.nifti.header(infiles[[1]])
    if (length(hdr$dim) == 2) {
        nvoxs <- hdr$dim[2]
    } else {
        nvoxs <- prod(hdr$dim[-length(hdr$dim)])
    }
    
    # Get mask
    if (is.null(mask)) {
        mask <- rep(TRUE, nvoxs)
    } else {
        mask <- read.mask(mask)
    }
    
    # Check match between functional and mask
    if (nvoxs != length(mask))
        stop("length mismatch between mask and functional: ", infiles[1])
    
    return(mask)
}

.get_nifti_nvols <- function(func_files1, func_files2=NULL, verbose=TRUE) {
    fun <- function(x) {
        hdr <- read.nifti.header(x)
        n <- length(hdr$dim)
        if (n == 4) {
            return(hdr$dim[4])
        } else if (n == 2) {
            return(hdr$dim[1])
        } else if (n == 1) {
            return(hdr$dim[1])
        } else {
            vstop("Input functional file '%s' must be 2 or 4 dimensions but is %i dimensional", x, n)
        }
    }
    
    progress <- ifelse(verbose, "text", "none")
    ntpts1 <- laply(func_files1, fun, .progress=progress)
    if (!is.null(func_files2)) {
        ntpts2 <- laply(func_files2, fun, .progress=progress)
        for (i in 1:length(ntpts1)) {
            if (ntpts1[i] != ntpts2[i]) {
                vstop("subject #%i does not have the same # of timepoints for the first and second functional datasets", i)
            }
        }
    }
    
    return(ntpts1)
}

.get_nifti_ndims <- function(func_files1, func_files2=NULL, verbose=TRUE) {
    fun <- function(x) {
        hdr <- read.nifti.header(x)
        length(hdr$dim)
    }
    
    progress <- ifelse(verbose, "text", "none")
    ndims1 <- laply(func_files1, fun, .progress=progress)
    ndims1[ndims1==1] <- 2
    if (!all(ndims1==ndims1[1]))
        stop("Not all dimensions in the first set of functionals are the same")
    if (ndims1[1] !=4 && ndims1[1] != 2)
        stop("Only allow 2 or 4 dimensional nifti files for first set of files")
    if (!is.null(func_files2)) {
        ndims2 <- laply(func_files2, fun, .progress=progress)
        ndims2[ndims2==1] <- 2
        if (!all(ndims2==ndims2[1]))
            stop("Not all dimensions in the second set of functionals are the same")
        if (ndims2[1] !=4 && ndims2[1] != 2)
            stop("Only allow 2 or 4 dimensional nifti files for second set of files")
    } else {
        ndims2 <- NULL
    }
    
    list(n1=ndims1[1], n2=ndims2[1])
}

.read_funcs <- function(infiles1, mask1, infiles2=NULL, mask2=NULL, 
                        verbose=TRUE, parallel=FALSE, shared=parallel, 
                        scale=TRUE, extra_checks=FALSE) 
{
    vcat(verbose, "...getting dimensions of functional data")
    ndims <- .get_nifti_ndims(infiles1, infiles2, verbose)
    
    # 1st Set of Functionals
    vcat(verbose, "Loading and masking functional data (Part 1)")
    ftype1 <- ifelse(ndims$n1==2, "nifti2d", "nifti4d")
    reader1 <- gen_big_reader(ftype1, type="double", shared=shared)
    funclist1 <- load_and_mask_func_data2(infiles1, reader1, mask=mask1, 
                                          verbose=verbose, scale=scale, 
                                          type="double", shared=shared)
    check1 <- check_func_data(infiles1[1], funclist1[1], extra=TRUE, 
                              verbose=verbose, parallel=FALSE)
    check2 <- check_func_data(infiles1[-1], funclist1[-1], extra=extra_checks, 
                              verbose=verbose, parallel=parallel)
    checks <- c(check1, check2)
    if (any(checks!=0)) {
        vcat(verbose, "Bad data for following files:")
        vcat(verbose, paste(infiles1[checks!=0], collapse="\n"))
        vstop("Quitting due to errors with 1st set of input functional data")
    }
    
    # 2nd Set of Functionals
    if ((!is.null(infiles2) && is.null(mask2)) || 
     (is.null(infiles2) && !is.null(mask2))) {
        stop("can't specify only infiles2 or only mask2")
    } else if (!is.null(infiles2) && !is.null(mask2)) {
        vcat(verbose, "Loading and masking functional data (Part 2)")
        ftype2 <- ifelse(ndims$n2==2, "nifti2d", "nifti4d")
        reader2 <- gen_big_reader(ftype2, type="double", shared=shared)
        funclist2 <- load_and_mask_func_data2(infiles2, reader2, mask=mask2,  
                                              verbose=verbose, scale=scale, 
                                              type="double", shared=shared)
        check1 <- check_func_data(infiles2[1], funclist2[1], extra=TRUE, 
                                  verbose=verbose, parallel=FALSE)
        check2 <- check_func_data(infiles2[-1], funclist2[-1], extra=extra_checks, 
                                  verbose=verbose, parallel=parallel)
        checks <- c(check1, check2)
        if (any(checks!=0)) {
            vcat(verbose, "Bad data for following files:")
            vcat(verbose, paste(infiles1[checks!=0], collapse="\n"))
            vstop("Quitting due to errors with 2nd set of input functional data")
        }
    } else {
        funclist2 <- NULL
    }
    
    return(list(funcs1=funclist1, funcs2=funclist2))
}
