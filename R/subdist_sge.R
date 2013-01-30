# subdist_sge_wrapper <- function()

# Returns design matrix
subdist.read_regressors <- function(design_fpath, verbose=TRUE)
{
    vcat(verbose, "Reading in regressors")
    
    if (!file.exists(design_fpath))
        vstop("Regressor file '%s' does not exist", design_fpath)
    
    design_mat <- as.matrix(read.table(design_fpath, header=TRUE))
    
    ## check rank deficiency
    k <- qlm_rank(as.big.matrix(design_mat))
    if (k < ncol(design_mat))
        vstop("design matrix '%s' is rank deficient", design_fpath)
    
    return(design_mat)
}

subdist.prepare_and_mask_funcs <- function(func_files, verbose=TRUE, ...)
{
    vcat(verbose, "Preparing and masking functionals")
    
    inlist <- load_funcs.prepare(func_files, verbose)
    inlist <- load_funcs.mask(inlist, verbose=verbose, ...)
    
    inlist$nsubs <- length(inlist1$files)
    inlist$nvoxs <- sum(inlist1$mask)
    
    return(inlist)
}

subdist.memory_limit <- function(opts, inlist1, inlist2=NULL, verbose=TRUE)
{
    vcat(verbose, "Determining memory limit/reqs")
    
    if (inlist1$nsubs != inlist2$nubs)
        stop("1st and 2nd set of functional data have different # of subjects")
    
    ntpts1 <- get_funclist_tpts(inlist1)
    
    if (!is.null(inlist2)) {
        ntpts2 <- get_funclist_tpts(inlist2)
        if (!all(ntpts1 == ntpts2)) {
            bad_subs <- which(ntpts1 != ntpts)
            bad_subs <- paste(bad_subs, collapse=", ")
            vstop("subjects %s do not have the same # of time-points for the 1st and 2nd functional datasets", bad_subs)
        }
    } else {
        nvoxs2 <- NULL
    }    
    
    nsubs <- inlist1$nsubs
    nvoxs1 <- inlist1$nvoxs
    nvoxs2 <- inlist2$nvoxs
    
    opts <- get_subdist_memlimit(opts, nsubs, nvoxs1, ntpts1, nvoxs2)
    
    return(opts)
}

subdist.check_funcs <- function(inlist, verbose=TRUE, extra_checks=FALSE, 
                                parallel=FALSE)
{
    vcat(verbose, "Checking functionals")
    
    # Detailed check on 1st subject and then a quick pass on the rest
    check1 <- check_func_data(inlist$files[1], inlist$funcs[1], extra=TRUE, 
                              verbose=verbose, parallel=FALSE)
    check2 <- check_func_data(inlist$files[-1], inlist$funcs[-1], 
                              extra=extra_checks, verbose=verbose, 
                              parallel=parallel)
    checks <- c(check1, check2)
    if (any(checks!=0)) {
        vcat(verbose, "Bad data for following files:")
        vcat(verbose, paste(inlist$files[checks!=0], collapse="\n"))
        vstop("Quitting due to errors with 1st set of input functional data")
    }
}

subdist.create_dists <- function(opts, outdir, inlist1, inlist2=NULL, ...)
{
    orig <- create_subdist(outdir, inlist1, inlist2, opts, ...)
    
    # Only keep the file paths
    ret <- list(sdist=orig$sdist.fname, gdist=orig$gdist.fname, bpath=orig$bpath)
    
    # Remove the actual distance matrices from memory
    rm(orig); invisible(gc(F,T))
    
    return(ret)
}

subdist_sge_wrapper <- function(funcfiles1, funcfiles2=NULL, design_mat=NULL, 
                                verbose=TRUE, extra_checks=FALSE, ...)
{
    if (!is.null(design_mat))
        design_mat <- subdist.read_regressors(design_mat, verbose)
    
    inlist1 <- subdist.prepare_and_mask_funcs(funcfiles1, verbose, ...)
    if (is.null(inlist2))
        inlist2 <- NULL
    else
        inlist2 <- subdist.prepare_and_mask_funcs(funcfiles2, verbose, ...)
    
    subdist.memory_limit(opts, inlist1, inlist2, verbose)
    
    subdist.check_funcs(inlist1, verbose, extra_checks=TRUE, 
                        parallel=FALSE)
    
}

# This script will prepare the relevant variables to run subdist_sge.runner
subdist_sge <- function(inlist1, list.dists, 
                              blocksize, superblocksize, 
                              inlist2==NULL, design_mat=NULL, 
                              verbose=1, parallel=FALSE)
{
    # Variables to indicate level of verbosity
    ## inform will show debugging information
    inform <- verbose==2
    verbosity <- verbose
    verbose <- as.logical(verbose)
    
    vcat(verbose, "Setup")
        
    # Super Blocks
    ## steps in which will go through the voxels or ROIs (related to inlist1)
    vcat(verbose, "...determine super-blocks")
    superblocks <- niftir.split.indices(1, nvoxs, by=superblocksize)
    
    
    
    preplist <- list(
        inform = inform, 
        verbosity = verbosity, 
        verbose = verbose, 
        dist_names = list(sdist=sdist.fname, gdist=gdist.fname, bpath=bpath), 
        nsubs = nsubs, 
        nvoxs = nvoxs, 
        design_mat = design_mat, 
        superblocks = superblocks
    )
}


