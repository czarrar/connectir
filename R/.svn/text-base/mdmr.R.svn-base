
#' Multi-Dimensional Multivariate Regression
#' 
#' Evaluates the degree to which a factor or predictor variable will explain 
#' the distances between particpants.
#' 
#' @author Zarrar Shehzad
#' @param G Matrix where each column represents one gower centered distance 
#'        matrix stored as a vector.
#' @param formula A formula specifying the model.
#' @param model A data frame in which the variables specified in the formula 
#'              will be found.
#' @param nperms Number of permutations (will add the original index)
#' @param factors2perm A vector indicating which factors to permute
#' @param superblocksize Number of voxels to process at a time
#' @param voxs A vector of voxel indices to examine
#' @param blocksize Number of permutations to process at a time
#' @param contr.unordered, contr.ordered Contrasts used for the design matrix
#' @param strata Vector of length nobs that contains indices to be exchanged
#' @param max.iter Maximum number of iterations for redoing a significantly 
#'                 correlated permutation
#' @param verbose 0, 1, 2 indicating nothing, moderate, a lot
#' @param parallel boolean
#' @param G.path Path to Gower's centered matrices if file-backed
#' @param fperms.path Path to Fperm's matrices
#' @peram save.fperms Boolean to save or not save the Fperms
#' @param type The data-type of the big matrices that will be created
#' @param shared Share intermediate and final big matrices
#' @param sge.info List containing njobs, nforks, nthreads, & ignore.proc.error
#' @param permute can be rhs, hat, or hat_with_covariate
#' @return list with modelinfo, pvals, fstats, and permutation indices
mdmr <- function(G, formula, model, 
                 nperms=4999, factors2perm=NULL, superblocksize=length(voxs), 
                 voxs=1:ncol(G), blocksize=nperms, 
                 contr.unordered="contr.sum", contr.ordered="contr.poly", 
                 strata=NULL, max.iter=100, 
                 verbose=1, parallel=FALSE, 
                 G.path=NULL, fperms.path=NULL, save.fperms=FALSE, 
                 type="double", shared=parallel, 
                 sge.info=NULL, permute="rhs")
{
    inform <- verbose==2
    verbose <- as.logical(verbose)
    progress <- ifelse(verbose, "text", "none")
    
    nvoxs <- length(voxs)
    nsubs <- sqrt(nrow(G))
    nperms <- nperms + 1    # account for original indices
    if ((superblocksize == nvoxs) && !is.null(sge.info))
        superblocksize <- round(superblocksize/2)
    superblocks <- niftir.split.indices(1, nvoxs, by=superblocksize)
    blocks <- niftir.split.indices(1, nperms, by=blocksize)
    blocks$nperms <- sapply(1:blocks$n, function(bi) {
                         blocks$ends[bi] - blocks$starts[bi] + 1
                     })
    
    # Checks
    if (!is.data.frame(model))
        stop("input model must be a data frame")
    if (!is.big.matrix(G) || !is.shared(G))
        stop("input Gower's matrix must be type big matrix and shared")
    if (!is.null(fperms.path) && !file.exists(fperms.path))
        stop("output path for Pseudo-F permutations does not exist")
    if (is.filebacked(G) && is.null(G.path))
        stop("G.path must be given when G is filebacked big matrix")
    if (!is.filebacked(G) && !is.null(G.path))
        warning("G.path will be ignored", immediate.=TRUE)
    if (!is.null(G.path) && !file.exists(G.path))
        stop("output path for distance matrices does not exist")
    if (nvoxs == 0)
        stop("number of voxels cannot be 0")
    if (nperms == 0)
        stop("number of permutations cannot be 0")    
    if (!is.null(sge.info)) {
        library(Rsge)
        ref.names <- c("njobs", "nforks", "nthreads", "ignore.proc.error")
        if (!all.equal(names(sge.info), ref.names))
            stop("sge.info must include njobs, nforks, nthreads, and ignore.proc.error")
        if (is.null(G.path))
            stop("must provide G.path for SGE MDMR job")
        if (save.fperms == T && is.null(fperms.path))
            stop("must provide fperms.path if saving for SGE MDMR job")
        sge.info$run <- TRUE
        sge.options(sge.user.options = sprintf("-S /bin/bash -pe mpi_smp %i", 
                                            sge.info$nthreads*sge.info$nforks))
    } else {
        sge.info <- list(run=FALSE)
    }
    
    # G filename
    G.desc <- describe(G)
    
    # Prepare model info
    # the list includes: rhs, qrhs, H2s, IHs, df.res, df.exp
    modelinfo <- mdmr_model(formula, model, contr.unordered, contr.ordered, 
                            factors2perm, verbose)
    factor.names <- names(attr(modelinfo$qrhs, "factors2perm"))
    nfactors <- length(factor.names)
    
    
    vcat(verbose, "Preparing data")
    
    # Permutations for each factor (fill in values later)
    vcat(verbose, "...preparing permutation matrices")
    list.perms <- mdmr_perms.gather_perms(modelinfo$rhs, modelinfo$qrhs, 
                                          nperms, strata, max.iter, 
                                          include.orig=TRUE, 
                                          verbose=verbose)
        
    # F-statistics (including permutations) for each factor
    # (will fill in values later)
    if (save.fperms) {
        vcat(verbose, "...preparing file-backed pesudo-F matrices")
        list.Fperms <- lapply(1:nfactors, function(fi) {
            if (is.null(fperms.path)) {
                bm <- big.matrix(nperms, nvoxs, type=type, shared=TRUE)
            } else {
                name <- factor.names[fi]
                bm <- big.matrix(nperms, nvoxs, type=type, shared=TRUE, 
                            backingpath=fperms.path, 
                            backingfile=sprintf("fperms_%s.bin", name), 
                            descriptorfile=sprintf("fperms_%s.desc", name))
            }
            return(bm)
        })
        list.fperms_desc <- lapply(list.Fperms, describe)
    } else {
        list.Fperms <- NULL
        list.fperms_desc <- NULL
    }
    
    vcat(verbose, "Computing MDMR across %i large blocks", 
         superblocks$n)
    vcat(verbose, "...with %i smaller blocks within each larger one",
         blocks$n)
    
    caller_for_superblocks <- function(si) {
        vcat(verbose, "large block %i", si)
        start.time <- Sys.time()
        
        if (sge.info$run == TRUE) {
            set_parallel_procs(sge.info$nforks, sge.info$nthreads, verbose, 
                               sge.info$ignore.proc.error)
        }
        
        firstVox <- superblocks$starts[si]
        lastVox <- superblocks$ends[si]
        sub_nvoxs <- lastVox - firstVox + 1
        
        # specify subset of G
        vcat(verbose, "...grabbing subset of Gower matrices")
        if (!is.null(G.path))
            G <- attach.big.matrix(G.desc, backingpath=G.path)
        tmp_Gs <- deepcopy(x=G, cols=firstVox:lastVox, type=type, shared=shared)
        if (!is.null(G.path))
            rm(G); invisible(gc(F,T))
        
        # prepare partial Fperms
        vcat(verbose, "...preparing partial pesudo-F matrices")
        list.partial_Fperms <- lapply(1:nfactors, function(fi) {
            big.matrix(nperms, sub_nvoxs, type=type, shared=TRUE)
        })
        
        # function to run mdmr for specific factor
        # with subset of permutations / voxels
        caller_for_mdmr_worker <- function(bi, include.orig=FALSE) {
            firstPerm <- blocks$starts[bi]
            lastPerm <- blocks$ends[bi]
            sub_nperms <- blocks$nperms[bi]
            include.orig <- ifelse(bi==1, TRUE, FALSE)
            
            # subset of partial Fmats
            vcat(verbose, "...extracting subset of partial f-stats")
            list.subset_Fperms <- lapply(list.partial_Fperms, function(partial_Fperms) {
                sub.big.matrix(partial_Fperms, firstRow=firstPerm, lastRow=lastPerm)
            })
            #list.subset_Fperms <- lapply(1:nfactors, function(fi) {
            #    big.matrix(sub_nperms, sub_nvoxs, type=type, shared=FALSE)
            #})
            
            # subset of permutation indices
            vcat(verbose, "...preparing subset of permutation indices")
            list.subset_perms <- lapply(1:nfactors, function(fi) {
                list.perms[[fi]][,firstPerm:lastPerm]
            })
            
            # calculate pseudo-F statistic
            # and return permutation indices
            mdmr.worker(modelinfo, tmp_Gs, list.subset_perms, 
                        list.subset_Fperms, permute=permute, 
                        type=type, shared=shared, 
                        verbose=verbose)
            
            # copy subset of partial Fmats
            #vcat(verbose, "...copying partial f-stats")
            #for (fi in 1:nfactors) {
            #    orig_Fperms <- list.subset_Fperms[[fi]]
            #    target_Fperms <- list.partial_Fperms[[fi]]
            #    bedeepcopy(x=orig_Fperms, y=target_Fperms, 
            #               y.rows=firstPerm:lastPerm)
            #}
            
            # remove subset Fmats and perm indices
            vcat(verbose, "...removing subset f-stats and permutation indices")
            rm(list.subset_Fperms, list.subset_perms); invisible(gc(F,T))
                       
            return(NULL)
        }
        
        # loop through blocks (permutations)
        vcat(verbose, "...looping through %i permutation block(s)", blocks$n)
        nulls <- llply(1:blocks$n, caller_for_mdmr_worker, 
                        .progress=progress, .inform=inform, 
                        .parallel=parallel)
        
        # calculate and save p-values
        vcat(verbose, "...calculating p-values")
        partial_Pmat <- mdmr.fstats_to_pvals(list.partial_Fperms)
        
        # save F-stats
        if (save.fperms) {
            vcat(verbose, "...saving F-statistics to disk")
            list.Fperms <- llply(1:nfactors, function(fi) {
                iFperms <- list.partial_Fperms[[fi]]
                
                if (!is.null(fperms.path)) {
                    Fperms <- attach.big.matrix(list.fperms_desc[[fi]], 
                                                backingpath=fperms.path)
                } else {
                    Fperms <- list.Fperms[[fi]]
                }
                oFperms <- sub.big.matrix(Fperms, 
                                firstCol=firstVox, lastCol=lastVox, 
                                backingpath=fperms.path)
                
                deepcopy(x=iFperms, y=oFperms)
                
                if (!is.null(fperms.path))
                    Fperms <- free.memory(Fperms, fperms.path)
                
                Fperms
            }, .progress=progress, .inform=inform, .parallel=parallel)
        }
        
        # remove temporary F-stats
        rm(list.partial_Fperms); invisible(gc(F,T))
        
        # time/duration
        end.time <- Sys.time()
        time.total <- as.numeric(end.time-start.time, units="mins")
        time.left <- time.total*(superblocks$n-si)
        vcat(verbose, "...took %.1f minutes (%.1f minutes left)\n", 
             time.total, time.left)
             
        return(partial_Pmat)
    }
    
    # Loop through super blocks (voxels)
    if (sge.info$run) {
        list.pvals <- sge.parLapply(1:superblocks$n, caller_for_superblocks, 
                                    debug=inform, trace=inform, 
                                    packages=c("connectir"), 
                                    function.savelist=ls(), 
                                    njobs=sge.info$njobs)
    } else {
        list.pvals <- lapply(1:superblocks$n, caller_for_superblocks)
    }
    
    # P-values for each factor
    vcat(verbose, "...compiling and saving p-values")
    if (is.null(fperms.path)) {
        Pmat <- big.matrix(nvoxs, nfactors, type=type, shared=FALSE)
    } else {
        Pmat <- big.matrix(nvoxs, nfactors, type=type, shared=TRUE, 
                           backingpath=fperms.path, 
                           backingfile="pvals.bin", 
                           descriptorfile="pvals.desc")
    }
    for (si in 1:superblocks$n) {
        firstVox <- superblocks$starts[si]
        lastVox <- superblocks$ends[si]
        oPmat <- list.pvals[[si]]
        Pmat[firstVox:lastVox,] <- oPmat
        #bedeepcopy(x=oPmat, y=Pmat, y.rows=firstVox:lastVox)
    }
    rm(list.pvals); invisible(gc(F,T))
    
    # Reload F-stats
    if (!is.null(fperms.path)) {
        vcat(verbose, "...reloading F-perms")
        rm(list.Fperms); invisible(gc(F,T))
        list.Fperms <- lapply(list.fperms_desc, function(desc) {
            attach.big.matrix(desc, backingpath=fperms.path)
        })
    }
    
    structure(list(
        modelinfo = modelinfo,   
        pvals = Pmat, 
        ns = list(nfactors=nfactors, nobs=nsubs, nperms=nperms, nvoxs=nvoxs), 
        fstats = list.Fperms, 
        fpath = fperms.path, 
        perms = list.perms
    ), class="mdmr")
}

save_mdmr <- function(obj, voxs, sdir, mdir, verbose=TRUE) {
    vcat(verbose, "Saving MDMR Results")
    
    mdmr.output <- mdir
    if (!file.exists(mdmr.output)) {
        vcat(verbose, "...creating MDMR output directory")
        dir.create(mdmr.output)
    }
    
    vcat(verbose, "...saving model info")
    save_mdmr.modelinfo(mdir, obj$modelinfo)
    
    vcat(verbose, "...saving formula")
    save_mdmr.formula(mdir, obj$modelinfo)
        
    vcat(verbose, "...saving factor names")
    save_mdmr.factors(mdir, obj$modelinfo)
    
    vcat(verbose, "...saving model")
    save_mdmr.model(mdir, obj$modelinfo)
    
    vcat(verbose, "...saving hat matrices")
    save_mdmr.hat_matrices(mdir, obj$modelinfo)
    
    vcat(verbose, "...saving permutations")
    factor.names <- names(attr(obj$modelinfo$qrhs, "factors2perm"))
    save_mdmr.permutations(mdir, obj$perms, factor.names)
    
    vcat(verbose, "...saving p-value and z-stat map(s)")
    save_mdmr.pvals_and_zstats(sdir, mdir, obj$pvals, voxs, factor.names)
    
    return(NULL)
}

save_mdmr.modelinfo <- function(mdir, modelinfo)
{
    ofile <- file.path(mdir, "modelinfo.rda")
    save(modelinfo, file=ofile)
}

save_mdmr.formula <- function(mdir, modelinfo)
{
    ofile <- file.path(mdir, "formula.txt")
    cat(as.character(modelinfo$formula)[-c(1:2)], file=ofile)
}

save_mdmr.factors <- function(mdir, modelinfo)
{
    ofile <- file.path(mdir, "factor_names.txt")
    
    all.factornames <- names(attr(modelinfo$qrhs, "factors.names"))
    perm.factornames <- names(attr(modelinfo$qrhs, "factors2perm"))
    is <- !(all.factornames %in% perm.factornames)
    cov.factornames <- all.factornames[is]
    
    totext <- sprintf("All Factors: %s\nPermuted Factors: %s\nCovariates: %s", 
                paste(all.factornames, collapse=", "), 
                paste(perm.factornames, collapse=", "), 
                paste(cov.factornames, collapse=", "))
    
    cat(totext, file=ofile)
}

save_mdmr.model <- function(mdir, modelinfo)
{
    ofile1 <- file.path(mdir, "model_raw.txt")
    ofile2 <- file.path(mdir, "model_evs.txt")
    
    write.csv(modelinfo$model, file=ofile1)
    write.csv(modelinfo$rhs, file=ofile2)
}

save_mdmr.hat_matrices <- function(mdir, modelinfo)
{
    factornames <- names(attr(modelinfo$qrhs, "factors2perm"))
    nfactors <- length(factornames)
    
    for (fi in 1:nfactors) {
        fname <- factornames[fi]
        ofile1 <- file.path(mdir, sprintf("H2s_%s.2D", fname))
        ofile2 <- file.path(mdir, sprintf("IHs_%s.2D", fname))
        write.table(modelinfo$H2s[[fi]], quote=F, row.names=F, col.names=F,  
                    file=ofile1)
        write.table(modelinfo$IHs[[fi]], quote=F, row.names=F, col.names=F,  
                    file=ofile2)
    }
}

save_mdmr.permutations <- function(mdir, list.perms, factor.names)
{
    nfactors <- length(list.perms)
    if (nfactors != length(factor.names)) {
        vstop("List of permutations must be length of %i and not %i", 
                nfactors, length(factor.names))
    }
    
    for (fi in 1:nfactors) {
        as.big.matrix(list.perms[[fi]], 
                      backingpath=mdir, 
                      descriptorfile=sprintf("perms_%s.desc", factor.names[fi]), 
                      backingfile=sprintf("perms_%s.bin", factor.names[fi]))
    }
}

save_mdmr.pvals_and_zstats <- function(sdir, mdir, Pvals, voxs, factor.names)
{
    maskfile <- file.path(sdir, "mask.nii.gz")
    mask <- read.mask(maskfile)
    hdr <- read.nifti.header(maskfile)

    nfactors <- ncol(Pvals)
    nvoxs <- nrow(Pvals)
    
    if (nfactors != length(factor.names))
        stop("Number of columns in Pvals must match length of factor.names")
    if (nvoxs != length(voxs))
        stop("Number of rows in Pvals must match length of voxs")
    
    for (i in 1:nfactors) {
        ofile1 <- file.path(mdir, sprintf("one_minus_pvals_%s.nii.gz", 
                                         factor.names[[i]]))
        ps <- vector("numeric", sum(mask))
        ps[voxs] <- 1 - Pvals[,i]
        write.nifti(ps, hdr, mask, outfile=ofile1, odt="float")
        
        ofile2 <- file.path(mdir, sprintf("zstats_%s.nii.gz", 
                                         factor.names[[i]]))
        zs <- vector("numeric", sum(mask))
        zs[voxs] <- qt(Pvals[,i], Inf, lower.tail=FALSE)
        write.nifti(zs, hdr, mask, outfile=ofile2, odt="float")
    }
}


#' INTERNAL: Runs MDMR
#' saving the p-value's and F-statistics
#' and returning a matrix of permuted observation indices
#' 
#' @author Zarrar Shehzad
#' @param modelinfo a list that includes rhs, QR of rhs, H2s, IHs, and dfs
#' @param mat.Gs matrix where each column is a Gower centered matrix
#' @param list.perms list of matrices of nobs x nperms
#' @param list.Fperms list of big matrices that will hold the F-statistics
#' @param verbose boolean
#' @param ... values given to \code{\link{mdmr_perms}}
#' @return list.Fperms
mdmr.worker <- function(modelinfo, mat.Gs, list.perms, list.Fperms, 
                        verbose=TRUE, ...)
{
    list.H2s <- mdmr_perms.gather_H2s(modelinfo$rhs, modelinfo$qrhs, 
                                      list.perms, verbose=verbose, ...)
    list.IHs <- mdmr_perms.gather_IHs(modelinfo$rhs, modelinfo$qrhs, 
                                      list.perms, verbose=verbose, ...)
    
    # fill in values for list.Fperms
    vcat(verbose, "...calculating Pseudo-F stats")
    .Call("mdmr_worker", mat.Gs, list.Fperms, list.H2s, list.IHs, 
            as.double(modelinfo$df.res), as.double(modelinfo$df.exp))
    
    return(list.Fperms)
}

#' INTERNAL: Calculate p-values from real and permuted F-stats
#' 
#' @author Zarrar Shehzad
#' @param list.Fperms list of big matrices that will hold the F-statistics
#' @param verbose boolean
#' @param parallel boolean
#' @return matrix of nvoxs x nfactors containing p-values
mdmr.fstats_to_pvals <- function(list.Fperms, verbose=TRUE, parallel=FALSE)
{
    nfactors <- length(list.Fperms)
    nvoxs <- ncol(list.Fperms[[1]])
    nperms <- nrow(list.Fperms[[1]])
    progress <- ifelse(verbose, "text", "none")
    
    for (fi in 1:nfactors) {
        if (ncol(list.Fperms[[fi]]) != nvoxs)
            vstop("# of voxels not all the same in list.Fperms")
        if (nrow(list.Fperms[[fi]]) != nperms)
            vstop("# of permutations not all the same in list.Fperms")
    }
    
    Pmat <- laply(1:nfactors, function(fi)
                    .Call("mdmr_fstats_to_pvals", list.Fperms[[fi]]), 
                    .progress=progress, .parallel=parallel)
    return(t(Pmat))
}
