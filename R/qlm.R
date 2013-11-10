# Quick Multiple Linear Regression using Big Matrices

# 1. Have a set of big matrices, for one voxel @ a time calculate a linear regression

# Given a design matrix => check it's rank; calculate the dd first + save
# want function that converts a bm => arma mat

setOldClass(c("qfit", "qcontrasts"))

formula_to_mat <- function(formula, data) {
    mf <- model.frame(formula=formula, data=data)
    X <- model.matrix(attr(mf, "terms"), data=mf)
    return(X)
}

formula_to_bmat <- function(formula, data) {
    mf <- model.frame(formula=formula, data=data)
    X <- model.matrix(attr(mf, "terms"), data=mf)
    return(as.big.matrix(X))
}

# Check if design matrix is rank deficient
setGeneric('qlm_rank', 
    function(X)
        standardGeneric('qlm_rank')
)

setMethod('qlm_rank',
    signature(X='big.matrix'),
    function(X) {
        .Call("big_qlm_rank", X, PACKAGE = "connectir")
    }
)


# Get diag( solve(t(X) %*% X) )
setGeneric('qlm_dd', 
    function(X)
        standardGeneric('qlm_dd')
)

setMethod('qlm_dd',
    signature(X='big.matrix'),
    function(X) {
        as.vector(.Call("big_qlm_dd", X, PACKAGE = "connectir"))
    }
)

# Solves for y ~ X
setGeneric('qlm_fit', 
    function(y, X, check.rank=TRUE, ...)
        standardGeneric('qlm_fit')
)

setMethod('qlm_fit',
    signature(y='big.matrix', X='big.matrix'),
    function(y, X, check.rank=FALSE, outputs=NULL, ...) {
        n = nrow(X); k = ncol(X); m = ncol(y)
        
        if (check.rank) {
            Xrank <- qlm_rank(X)
            if (Xrank < k)
                stop("design matrix is rank deficient.")
        }
        
        if (is.null(outputs)) {
            outputs <- list(
                coefficients = big.matrix(k, m, type="double", ...), 
                residuals = big.matrix(nrow(y), m, type="double", ...), 
                mse = big.matrix(m, 1, type="double", ...)
            )
        } else {
            if (!all(c("coefficients", "residuals", "mse") %in% outputs))
                stop("outputs must be list of 'coefficients', 'residuals', and 'mse'.")
            if (nrow(outputs$coef) != n && ncol(outputs$coef) != m)
                stop("size mismatch for output coefficients")
            if (nrow(outputs$res) != n && ncol(outputs$res) != m)
                stop("size mismatch for output residuals")
            if (nrow(outputs$mse) != m && ncol(outputs$mse) != 1)
                stop("size mismatch for output mse")
        }
        
        .Call("big_qlm_fit", y, X, outputs$coef, outputs$res, outputs$mse, 
                as.double(n), as.double(k), as.double(m), PACKAGE = "connectir")
        
        outputs$n <- n; outputs$k <- k; outputs$m <- m
        outputs$call <- match.call()
        class(outputs) <- "qfit"
        
        return(outputs)
    }
)

# R^2 (after qlm_fit)
setGeneric('qlm_rsquared', 
    function(y, inputs, ...)
        standardGeneric('qlm_rsquared')
)

setMethod('qlm_rsquared',
    signature(y='big.matrix', inputs='qfit'),
    function(y, inputs, adjust=TRUE, output=NULL, ...) {
        if (is.null(output)) {
            output <- big.matrix(inputs$m, 1, type="double", ...)
        } else if (nrow(output) != inputs$n && ncol(output) != inputs$m) {
            stop("size mismatch for output matrix")
        }
        .Call("big_qlm_rsquared", y, inputs$coef, inputs$res, inputs$mse, 
                as.double(inputs$n), as.double(inputs$k), 
                output, as.integer(adjust), 
                PACKAGE = "connectir")
        
        return(output)
    }
)

# Residuals
setGeneric('qlm_residuals', 
    function(y, X, check.rank=FALSE, residuals=NULL, add.mean=FALSE, ...)
        standardGeneric('qlm_residuals')
)

setMethod('qlm_residuals',
    signature(y='big.matrix', X='big.matrix'),
    function(y, X, check.rank=FALSE, residuals=NULL, add.mean=FALSE, ...) {
        n = nrow(X); k = ncol(X); m = ncol(y)
        
        if (check.rank) {
            Xrank <- qlm_rank(X)
            if (Xrank < k)
                stop("design matrix is rank deficient.")
        }
        
        if (is.null(residuals)) {
            residuals = big.matrix(nrow(y), ncol(y), type="double", ...)
        } else if (nrow(residuals) != n && ncol(residuals) != m) {
            stop("size mismatch for output residuals")
        }
        
        .Call("big_qlm_residuals", y, X, residuals, as.integer(add.mean), 
				PACKAGE = "connectir")
        
        return(residuals)
    }
)

# Gets tvals for given contrasts from qlm_fit results
setGeneric('qlm_contrasts', 
    function(inputs, contrasts, dd, ...)
        standardGeneric('qlm_contrasts')
)

setMethod('qlm_contrasts',
    signature(inputs='qfit', contrasts='matrix', dd='vector'),
    function(inputs, contrasts, dd, outputs=NULL, ...) {
        if (ncol(contrasts) != inputs$k)
            stop("column mismatch for contrasts")
        if (length(dd) != inputs$k)
            stop("length mismatch for dd")
        n <- nrow(contrasts)
        
        if (is.null(outputs)) {
            outputs <- list(
                coefficients = big.matrix(n, inputs$m, type="double", ...), 
                standard_errors = big.matrix(n, inputs$m, type="double", ...), 
                tvalues = big.matrix(n, inputs$m, type="double", ...)
            )
        } else {
            if (!all(c("coefficients", "standard_errors", "tvalues") %in% outputs))
                stop("outputs must be list of 'coefficients', 'standard_errors', and 'tvalues'")
            if (nrow(outputs$coef) != n && ncol(outputs$coef) != inputs$m)
                stop("size mismatch for output coefficients")
            if (nrow(outputs$st) != n && ncol(outputs$st) != inputs$m)
                stop("size mismatch for output standard errors")
            if (nrow(outputs$tval) != n && ncol(outputs$tval) != inputs$m)
                stop("size mismatch for output t-values")
        }
        
        .Call("big_qlm_contrasts", inputs$coef, inputs$mse, contrasts, dd, outputs$coef, 
                outputs$st, outputs$tval, as.double(inputs$m), PACKAGE = "connectir")
        
        outputs$n <- n
        outputs$call <- match.call()
        class(outputs) <- "qcontrasts"
        
        return(outputs)
    }
)

# Checks design matrix
# and if rank deficient removes offending columns
check_design_matrix <- function(X, to.quit=FALSE, verbose=TRUE) {
    qrhs <- qr(X)
    if (qrhs$rank < ncol(X)) {
        if (to.quit)
            stop("model is rank deficient")
        
        vcat(verbose, " warning: model is rank deficient")
        bad_cols <- paste(qrhs$pivot[-c(1:qrhs$rank)], collapse=", ")
        bad_colnames <- paste(colnames(X)[bad_cols], collapse=", ")
        vcat(verbose, " will drop cols %s (%s) from design matrix", 
                        bad_cols, bad_colnames)
        
        X <- X[, X$pivot, drop = FALSE]
        X <- X[, 1:X$rank, drop = FALSE]
    }
    
    return(X)
}

vox_glm <- function(funclist, evs, cons, blocksize, outmats, bp=NULL, 
                    verbose=TRUE, parallel=FALSE, shared=parallel, 
                    ztransform=FALSE)
{
    if (!is.matrix(cons))
        stop("cons must be a matrix")
    vcat(verbose, "...setup")
    nsubs <- length(funclist)
    nvoxs <- ncol(funclist[[1]])
    voxs <- 1:nvoxs
    ncons <- nrow(cons)
    blocks <- niftir.split.indices(1, nvoxs, by=blocksize)
    progress <- ifelse(verbose, "text", "none")
    k <- qlm_rank(evs)
    if (k < ncol(evs))
        stop("EV matrix is rank deficient")
    dd <- qlm_dd(evs)
    if (!is.shared(outmats[[1]]))
        stop("outmats must be shared")
    if (is.filebacked(outmats[[1]]) && is.null(bp))
        stop("backingpath (bp) must be given if outmats are file-backed")
    if (length(outmats) != ncons)
        stop("cons doesn't match with outmats")
    
    gfun <- function(bi, outmats) {
        vcat(verbose, "...block %i", bi)
        
        first <- blocks$starts[bi]; last <- blocks$ends[bi]
        sub_voxs <- first:last
        sub_nvoxs <- length(sub_voxs)
        
        tmp_outmats <- lapply(1:ncons, function(x) {
                            big.matrix(nvoxs, sub_nvoxs, init=0, type="double", shared=shared)
                        })
        
        # correlate
        vcat(verbose, "....correlate")
        subs.cormaps <- vbca_batch2(funclist, c(first, last), 
                                    ztransform=ztransform, 
                                    type="double", shared=shared)
        
        # loop through each seed
        vcat(verbose, '....glm')
        sfun <- function(si) {
            # combine and transpose correlation maps for 1 seed
            seedCorMaps <- big.matrix(nsubs, nvoxs-1, type="double", shared=FALSE)
            .Call("subdist_combine_and_trans_submaps", subs.cormaps, as.double(si), 
                  as.double(voxs[-sub_voxs[si]]), seedCorMaps, PACKAGE="connectir")
                        
            # glm
            fit <- qlm_fit(seedCorMaps, evs, shared=FALSE)
            res <- qlm_contrasts(fit, cons, dd, shared=FALSE)
            
            for (ci in 1:ncons)
                tmp_outmats[[ci]][-voxs[sub_voxs[si]],si] <- res$tvalues[ci,]
            
            rm(seedCorMaps, fit, res); gc(FALSE)
            
            return(NULL)
        }
        
        llply(1:sub_nvoxs, sfun, .parallel=parallel, .progress=progress, .inform=T)
        
        # copy over temp data
        vcat(verbose, '....save')
        for (ci in 1:ncons) {
            sub_outmat <- sub.big.matrix(outmats[[ci]], firstCol=first, lastCol=last, backingpath=bp)
            deepcopy(x=tmp_outmats[[ci]], y=sub_outmat)
            if (!is.null(bp)) {
                flush(sub_outmat); flush(outmats[[ci]])
                outmats[[ci]] <- free.memory(outmats[[ci]], bp)
            }
        }
        
        rm(subs.cormaps, tmp_outmats); gc(FALSE)
        
        return(outmats)
    }
    
    vcat(verbose, "...%i blocks", blocks$n)
    for (bi in 1:blocks$n)
        outmats <- gfun(bi, outmats)
    
    invisible(outmats)
}

vox_glm3 <- function(funclist1, evs, cons, blocksize, 
                      funclist2=funclist1, 
                      outmats=NULL, bp=NULL, 
                      verbose=TRUE, parallel=FALSE, shared=parallel, 
                      ztransform=FALSE)
{
    vcat(verbose, "...checks")
    if (is.null(funclist2)) funclist2 <- funclist1
    if (!is.matrix(cons))
        stop("cons must be a matrix")
    if (!is.null(outmats) && !is.shared(outmats[[1]]))
        stop("outmats must be shared")
    if (!is.null(outmats) && (is.filebacked(outmats[[1]]) && is.null(bp)))
        stop("backingpath (bp) must be given if outmats are file-backed")
    if (length(funclist1) != length(funclist2))
        stop("length mismatch between first and second set of functionals")
    
    vcat(verbose, "...setup")
    nsubs <- length(funclist1)
    nvoxs1 <- ncol(funclist1[[1]])
    nvoxs2 <- ncol(funclist2[[1]])
    voxs1 <- 1:nvoxs1
    voxs2 <- 1:nvoxs2
    ncons <- nrow(cons)
    blocks <- niftir.split.indices(1, nvoxs1, by=blocksize)
    progress <- ifelse(verbose, "text", "none")
    k <- qlm_rank(evs)
    if (k < ncol(evs))
        stop("EV matrix is rank deficient")
    dd <- qlm_dd(evs)

    if (is.null(outmats)) {
        vcat(verbose, "...creating output matrices")
        outmats <- lapply(1:ncons, function(i) {
            bfile <- sprintf("tvals_%02i.bin", i)
            dfile <- sprintf("tvals_%02i.desc", i)
            big.matrix(nvoxs2, nvoxs1, type="double", 
                       backingfile=bfile, descriptorfile=dfile, 
                       backingpath=bp)
        })
    } else {
        if (length(outmats) != ncons)
            stop("cons doesn't match with outmats")
        if (nrow(outmats[[1]]) != nvoxs2)
            stop("incorrect # of rows for outmats")
        if (ncol(outmats[[1]]) != nvoxs1)
            stop("incorrect # of cols for outmats")
    }
    
    gfun <- function(bi, outmats) {
        vcat(verbose, "...block %i", bi)
        
        first <- blocks$starts[bi]; last <- blocks$ends[bi]
        sub_voxs <- first:last
        sub_nvoxs <- length(sub_voxs)
        
        tmp_outmats <- lapply(1:ncons, function(x) {
                            big.matrix(nvoxs2, sub_nvoxs, init=0, 
                                        type="double", shared=shared)
                        })
        
        # correlate
        vcat(verbose, "....correlate")
        subs.cormaps <- vbca_batch3(funclist1, c(first, last), funclist2, 
                                    ztransform=ztransform, 
                                    type="double", shared=shared)
        
        # loop through each seed
        vcat(verbose, '....glm')
        sfun <- function(si) {
            # combine and transpose correlation maps for 1 seed
            seedCorMaps <- big.matrix(nsubs, nvoxs2, type="double", shared=FALSE)
            .Call("subdist_combine_and_trans_submaps", subs.cormaps, as.double(si), 
                  as.double(voxs2), seedCorMaps, PACKAGE="connectir")
            
            # glm
            fit <- qlm_fit(seedCorMaps, evs, shared=FALSE)
            res <- qlm_contrasts(fit, cons, dd, shared=FALSE)
            
            for (ci in 1:ncons)
                tmp_outmats[[ci]][,si] <- res$tvalues[ci,]
            
            rm(seedCorMaps, fit, res); gc(FALSE)
            
            return(NULL)
        }
        
        llply(1:sub_nvoxs, sfun, .parallel=parallel, .progress=progress, .inform=T)
        
        # copy over temp data
        vcat(verbose, '....save')
        for (ci in 1:ncons) {
            sub_outmat <- sub.big.matrix(outmats[[ci]], firstCol=first, lastCol=last, 
                                         backingpath=bp)
            deepcopy(x=tmp_outmats[[ci]], y=sub_outmat)
            if (!is.null(bp)) {
                flush(sub_outmat); flush(outmats[[ci]])
                outmats[[ci]] <- free.memory(outmats[[ci]], bp)
            }
        }
        
        rm(subs.cormaps, tmp_outmats); gc(FALSE)
        
        return(outmats)
    }
    
    vcat(verbose, "...%i blocks", blocks$n)
    for (bi in 1:blocks$n)
        outmats <- gfun(bi, outmats)
    
    invisible(outmats)
}

# glm_summarize
# get average at each voxel
# get max at each voxel
# get number of voxels that are significant
# get number of voxels that are significant after fdr correction
# get standard deviation
glm_summarize <- function(tmats, tdir, df, pthr=0.05, 
                          outdir=NULL, outhdr=NULL, outmask=NULL, 
                          to_return=FALSE, 
                          memlimit=4, 
                          verbose=TRUE, parallel=FALSE, shared=parallel)
{
    vcat(verbose, "Summarizing GLM Results")
    measures <- c("mean_tstat", "sd_tstat", "max_tstat", "min_tstat", 
                  sprintf("percent_signif_%f_tstat", pthr))
    nmeasures <- length(measures)
    thresh <- qt(pthr/2, df, lower.tail=F)
    
    if (is.null(outdir) && !to_return)
        stop("must specify at least one of outdir or to_return")
    if (!is.null(outdir) && !file.exists(outdir))
        vstop("output directory '%s' doesn't exist", outdir)
    test <- c(!is.null(outdir), !is.null(outhdr), !is.null(outmask))
    if (any(test) && !all(test))
        stop("must specify none or all: outdir, outhdr, outmask")
    if (all(test))
        to_save <- TRUE
    else
        to_save <- FALSE
    
    vcat(verbose, "...re-reading file-backed data")
    tmp_tmats <- lapply(tmats, function(x) {
        fname <- paste(rmext(describe(x)@description$filename), ".desc", sep="")
        attach.big.matrix(file.path(tdir, fname))
    })
    n <- length(tmats)
    nr <- nrow(tmp_tmats[[1]])
    nc <- ncol(tmp_tmats[[1]])
    
    vcat(verbose, "...determining memory demands given limit of %.2f GB", memlimit)
    nforks <- getDoParWorkers()
    mem_out <- n2gb(n*nc*nmeasures)
    mem_col <- n2gb(nr*nforks)
    f <- function(x) memlimit - mem_out - mem_col*x
    m <- f(nc)
    if (m > 0) {
        blocksize <- nc
    } else {
        vcat(verbose, paste("...if you wanted to hold everything in memory",
                                 "you would need at least %.2f GB of RAM"), 
                           memlimit - m)
        vcat(verbose, "...autosetting blocksize")
        x <- tryCatch(floor(uniroot(f, c(1,nc))$root), error=function(ex) NA)
        if (is.na(x))
            stop("You don't have enough RAM")
        blocksize <- x
    }
    vcat(verbose, "...setting block size to %i (out of %i columns)", 
         blocksize, nc)
    if (blocksize < 1)
        stop("block size is less than 1")
    blocks <- niftir.split.indices(1, nc, by=blocksize)
    
    vcat(verbose, "...creating output vectors/matrices")
    outmats <- lapply(1:n, function(i) {
        tmp <- matrix(NA, nc, nmeasures)
        colnames(tmp) <- measures
        tmp
    })
    
    vcat(verbose, "...summarizing")
    sfun <- function(ci) {
        inds <- sblocks$starts[ci]:sblocks$ends[ci]
        # 1. Average
        outmat[inds,1] <- colmean(tmat, inds)
        # 2. Standard Deviation
        outmat[inds,2] <- colsd(tmat, inds)
        # 3. Max
        outmat[inds,3] <- colmax(tmat, inds)
        # 4. Min
        outmat[inds,4] <- colmin(tmat, inds)
        # 5. Percent of significant voxels
        outmat[inds,5] <- colMeans(abs(tmat[,inds])>thresh)
    }
    for (ti in 1:n) {
        vcat(verbose, "...tmat #%i with %i blocks", ti, blocks$n)
        tmat <- tmats[[ti]]
        outmat <- outmats[[ti]]
        for (bi in 1:blocks$n) {
            sblocks <- niftir.split.indices(blocks$starts[bi], blocks$ends[bi], 
                                            by=nforks)
            laply(1:sblocks$n, sfun, .inform=verbose, .parallel=parallel)
            tmat <- free.memory(tmats[[ti]], tdir)
        }
    }
    
    if (to_save) {
        vcat(verbose, "...saving")
        for (ti in 1:n) {
            outmat <- outmats[[ti]]
            for (mi in 1:nmeasures) {
                outfn <- file.path(outdir, sprintf("%s_%02i.nii.gz", measures[mi], ti))
                write.nifti(outmat[,mi], outhdr, outmask, outfile=outfn, odt="float")
            }
        }
    }
    
    rm(tmp_tmats); gc(FALSE, TRUE)
    
    if (to_return) {
        vcat(verbose, "...returning")
        return(outmats)
    }
}

