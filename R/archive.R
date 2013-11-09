.archive.load_and_mask_func_data <- function(fnames, mask, check=TRUE, type=NULL, 
                                             verbose=FALSE, ...) 
{
    if (is.character(mask))
        mask <- read.mask(mask)
    if (is.character(fnames) && !is.vector(fnames))
        fnames <- list(fnames)
    
    if (verbose)
        progress="text"
    else
        progress="none"
    
    if (verbose)
        cat("reading in data\n")
    dat.list <- llply(fnames, function(f) {
        if (is.null(type))
            x <- read.big.nifti4d(f, ...)
        else
            x <- read.big.nifti4d(f, type=type, ...)
        x <- do.mask(x, mask)
        y <- scale(x, to.copy=T, ...) # fix with update of deepcopy!
        rm(x)
        gc(FALSE)
        y
    }, .progress=progress)
    
    if (verbose)
        cat("checking data\n")
    if (check) {
        nc <- ncol(dat.list[[1]])
        l_ply(1:length(dat.list), function(i) {
            dat <- dat.list[[i]]
            fname <- fnames[[i]]
            # Everything must have same # of voxels
            if (nc != ncol(dat))
                stop(sprintf("%s must have the same # of nodes (%i vs %i) as other datasets", 
                                fname, nc, ncol(dat)))
            
            # There can't be any NaNs
            col.nas <- colna(dat.list[[i]])>0
            if (any(col.nas)) {
                w <- paste(which(col.nas), collapse=", ")
                stop(sprintf("%s has NaNs in nodes %s", fname, w))
            }
        }, .progress=progress)
    }
    
    gc(FALSE)
    
    if (length(fnames) == 1)
        return(dat.list[[1]])
    else
        return(dat.list)
}
