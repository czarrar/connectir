# subdist_sge_wrapper <- function()

# Returns design matrix
subdist.read_regressors <- function(fpath, formula=NULL, verbose=TRUE)
{
    vcat(verbose, "Reading in regressors")
    
    if (!file.exists(fpath))
        vstop("regressor file '%s' does not exist", fpath)
    
    model <- read.table(fpath, header=TRUE)
    
    if (any(is.na(model)) || any(is.null(model)))
        vstop("There is a missing element in regressor file '%s'", fpath)
    
    # Convert to design matrix if input is data frame
    if (!is.null(formula)) {
        rhs.frame <- model.frame(formula, model, drop.unused.levels = TRUE)
        rhs <- model.matrix(formula, rhs.frame)
    } else {
        rhs <- model
    }
    
    rhs <- as.matrix(rhs[,])
    
    ## check rank deficiency
    k <- qlm_rank(as.big.matrix(rhs))
    if (k < ncol(rhs))
        vstop("regressor design matrix '%s' is rank deficient", fpath)
    
    return(rhs)
}

subdist.prepare_funcs <- function(func_files, verbose=TRUE, ...) {
    vcat(verbose, "Preparing functionals")
    inlist <- load_funcs.prepare(func_files, verbose, ...)
    inlist$nsubs <- length(inlist$files)
    return(inlist)
}

subdist.prepare_mask <- function(inlist, verbose=TRUE, ...) {
    inlist <- load_funcs.mask(inlist, verbose=verbose, ...)
    inlist$nvoxs <- sum(inlist$mask)
    return(inlist)
}

subdist.prepare_and_mask_funcs <- function(func_files, verbose=TRUE, type="double", ...)
{
    vcat(verbose, "Preparing and masking functionals")
    
    inlist <- load_funcs.prepare(func_files, verbose, type)
    inlist <- load_funcs.mask(inlist, verbose=verbose, ...)
    
    inlist$nsubs <- length(inlist$files)
    inlist$nvoxs <- sum(inlist$mask)
    
    return(inlist)
}

# TODO: redo memory 
subdist.memory_limit <- function(memlimit, blocksize, superblocksize, 
                                 inlist1, inlist2=NULL, verbose=TRUE, 
                                 nforks=NULL)
{
    vcat(verbose, "Memory?")
    
    ntpts1 <- get_funclist_tpts(inlist1)
    
    if (!is.null(inlist2)) {
        if (inlist1$nsubs != inlist2$nubs)
            stop("1st and 2nd set of functional data have different # of subjects")
        
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
    
    opts <- list(memlimit=memlimit, verbose=verbose, 
                 blocksize=blocksize, superblocksize=superblocksize)
    opts <- get_subdist_memlimit(opts, nsubs, nvoxs1, ntpts1, nvoxs2, nforks)
    
    list(blocksize=opts$blocksize, superblocksize=opts$superblocksize)
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

subdist.create_dists <- function(opts, outdir, inlist1, inlist2=NULL, ret.orig=F, ...)
{
    orig <- create_subdist(outdir, inlist1, inlist2, opts, ...)
    
    if (ret.orig) {
        return(orig)
    } else {
        # Only keep the file paths
        ret <- list(sdist=orig$sdist.fname, gdist=orig$gdist.fname, bpath=orig$bpath)
        # Remove the actual distance matrices from memory
        rm(orig); invisible(gc(F,T))
        return(ret)
    }
}

subdist_wrapper <- function(funcfiles1, funcfiles2=NULL, design_mat=NULL, 
                            verbose=TRUE, parallel=FALSE, 
                            extra_checks=FALSE, 
                            # sge.info, memlimit, blocksize, superblocksize, nforks
                            ...)
{
    ###
    # Setup
    ###
    
    # 1. Covariates of Non-Intererest
    if (!is.null(design_mat))
        design_mat <- subdist.read_regressors(design_mat, verbose)
    
    # 2. Prepare functional files and mask
    inlist1 <- subdist.prepare_and_mask_funcs(funcfiles1, verbose, ...)
    if (is.null(funcfiles2))
        inlist2 <- NULL
    else
        inlist2 <- subdist.prepare_and_mask_funcs(funcfiles2, verbose, ...)
    
    # 3. Check memory limits
    ###  and Set block sizes
    l <- subdist.memory_limit(memlimit, blocksize, superblocksize, 
                              inlist1, inlist2, verbose, nforks)
    blocksize <- l$blocksize
    superblocksize <- l$superblocksize
    superblocks <- niftir.split.indices(1, inlist1$nvoxs, by=superblocksize)
    
    # 4. Check input functionals
    subdist.check_funcs(inlist1, verbose, extra_checks=TRUE, 
                        parallel=parallel)
    if (!is.null(inlist2)) {
        subdist.check_funcs(inlist2, verbose, extra_checks=extra_checks, 
                            parallel=parallel)
    }
    
    # 5. Setup distances
    list.dists <- subdist.create_dists(opts, outdir, inlist1, inlist2, ...)
    
    # 6. Scale the time-series?
    ## no scaling if connectivity is computed via an inverse covariance matrix
    glasso <- list(...)$glasso
    scale <- ifelse(is.null(glasso), FALSE, !glasso)
    
    
    ###
    # When not using sge
    ###
    
    if (is.null(sge.info)) {
        # Read and scale input functionals
        inlist1 <- load_funcs.read_and_scale(inlist1, verbose, to.copy=FALSE, 
                                             parallel=parallel, scale=scale, 
                                             ...)
        if (!is.null(inlist2)) {
            inlist2 <- load_funcs.read_and_scale(inlist2, verbose, 
                                                 to.copy=FALSE, 
                                                 parallel=parallel, 
                                                 scale=scale, ...)
        }
        
        # Call regular subdist wrapper
        # TODO
    }
    
    # Call subdist sge wrapper
    
    subdist_sge(inlist1, list.dists, blocksize, superblocksize, inlist2, 
                design_mat, verbose*1, parallel, ...)
}

