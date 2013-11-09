# BELOW: CODE DIRECTLY TAKEN FROM vegan
#permuted.index <- function (n, strata) 
#{
#    if (missing(strata) || is.null(strata)) 
#        out <- sample(n, n)
#    else {
#        out <- 1:n
#        inds <- names(table(strata))
#        for (is in inds) {
#            gr <- out[strata == is]
#            if (length(gr) > 1) 
#                out[gr] <- sample(gr, length(gr))
#        }
#    }
#    out
#}
# ABOVE: CODE DIRECTLY TAKEN FROM vegan

# NOTE: assume that y-intercept exists...need to deal with case when don't have intercept
mdmr.prepare.model <- function(formula, model, 
                               contr.unordered="contr.sum", 
                               contr.ordered="contr.poly", 
                               factors2perm=NULL, 
                               verbose=FALSE) 
{
    vcat(verbose, "...preparing model and stuff")
    TOL <- 1e-07

    n <- nrow(model)
    Terms <- terms(formula, data = model)
    formula[[2]] <- NULL
    
    # Get right hand matrix
    rhs.frame <- model.frame(formula, model, drop.unused.levels = TRUE)
    if (nrow(rhs.frame) != nrow(model))
        vstop("One of your factors has an empty element")
    op.c <- options()$contrasts
    options(contrasts = c(contr.unordered, contr.ordered))
    rhs <- model.matrix(formula, rhs.frame)
    options(contrasts = op.c)
    grps <- attr(rhs, "assign")
    qrhs <- qr(rhs)
    
    # deal with rank deficiency
    if (qrhs$rank < ncol(rhs)) {
        bad_cols <- paste(qrhs$pivot[-c(1:qrhs$rank)], collapse=", ")
        bad_colnames <- paste(colnames(rhs)[bad_cols], collapse=", ")
        vcat(verbose, " warning: model is rank deficient")
        vcat(verbose, " will drop cols %s (%s) from design matrix", 
                        bad_cols, bad_colnames)
    }
    rhs <- rhs[, qrhs$pivot, drop = FALSE]
    rhs <- rhs[, 1:qrhs$rank, drop = FALSE]
    
    # Get different groups of variables
    grps <- grps[qrhs$pivot][1:qrhs$rank]
    u.grps <- unique(grps)
    nterms <- length(u.grps) - 1
    factor.names <- attr(attr(rhs.frame, "terms"), "term.labels")[u.grps]
    
    # Get hat matrices
    H2s <- lapply(2:length(u.grps), function(j) {
        Xj <- rhs[, grps %in% u.grps[1:j]]
        qrX <- qr(Xj, tol = TOL)
        Q <- qr.Q(qrX)
        tcrossprod(Q[, 1:qrX$rank])
    })
    H <- H2s[[nterms]]  # complete model
    IH <- diag(n) - H
    if (length(H2s) > 1) 
        for (i in length(H2s):2) 
            H2s[[i]] <- H2s[[i]] - H2s[[i - 1]]
    
    # Get degrees of freedom
    df.Res <- n - qrhs$rank
    df.Exp <- sapply(u.grps[-1], function(i) sum(grps == i))
    
    # Factors to permute (i.e., assess significance)
    if (is.null(factors2perm)) {
        factors2perm <- 1:length(factor.names)
    } else if (is.character(factors2perm)) {
        factors2perm <- sapply(factors2perm, function(x) 
                                which(factor.names==x))
    }
    factor2perm.names <- factor.names[factors2perm]
    ## refine degrees of freedom
    df.Exp <- sapply(factors2perm, function(i) df.Exp[i])
    
    return(list(
        rhs = rhs,
        grps = grps,
        u.grps = u.grps, 
        factor.names = factor.names, 
        factors2perm = factors2perm, 
        factor2perm.names = factor2perm.names, 
        H = H,
        H2s = H2s,
        IH = IH,
        df.Res = df.Res,
        df.Exp = df.Exp,
        nobs = n,
        nfactors = length(factor.names), 
        nfactors2perm = length(factor2perm.names), 
        verbose = verbose, 
        progress = ifelse(verbose, "text", "none")
    ))
}

mdmr.prepare.permutations <- function(modelinfo, nperms, strata, max.iter) {
    vcat(modelinfo$verbose, "...preparing permutations")

    # get matrix containing different permutations of observations
    p <- sapply(1:nperms, function(x) permuted.index(modelinfo$nobs, strata))
    
    # remove permutations that are significantly correlated with your model
    if (max.iter > 0) {
        # r value that is significant (p < 0.05)
        rthresh <- tanh(1.96/sqrt(modelinfo$nobs-3))
        which.cor <- 1:nperms
        for (iter in 1:max.iter) {
            # get correlations between real model and permuted ones
            which.cor <- unique(unlist(lapply(modelinfo$factors2perm, function(i) {
                H <- modelinfo$H2s[[i]]
                tmp <- sapply(which.cor, function(j) cor(H[,1], H[p[,j],1]))
                which(abs(tmp)>rthresh)
            })))
            
            # if no significant correlations, break
            if (length(which.cor)==0)
                break
            
            # get new permutations to replace correlated ones
            p[,which.cor] <- sapply(1:length(which.cor), function(x) 
                                    permuted.index(modelinfo$nobs, strata))
        }
    }
    
    if (modelinfo$nobs < 2^31) {
        dims <- dim(p)
        p <- as.integer(p)
        dim(p) <- dims
    }
    
    return(p)
}

mdmr.prepare.permH2s <- function(modelinfo, p, firstRow, lastRow, ...) {
    if (firstRow > lastRow)
        stop("first row can't be greater than the last row")
    prows <- firstRow:lastRow
    nperms <- length(prows)
    if (nperms > ncol(p))
        stop("more rows requested than permutations")
    
    nfactors <- length(modelinfo$factors2perm)
    
    H2mats <- lapply(1:nfactors, function(ii) {
        vcat(modelinfo$verbose, "...preparing permuted model matrices for factor #%i", 
             ii)
        
        i <- modelinfo$factors2perm[ii]
        bigmat <- big.matrix(modelinfo$nobs^2, nperms+1, ...)
        bigmat[,1] <- as.vector(modelinfo$H2s[[i]])
        
        l_ply(1:nperms, function(kk) {
            k <- prows[kk]
            prm <- p[,k]
            bigmat[,kk+1] <- as.vector(modelinfo$H2s[[i]][prm,prm])
        }, .progress=modelinfo$progress)
        
        return(bigmat)
    })
    
    return(H2mats)
}

mdmr.prepare.permIH <- function(modelinfo, p, firstRow, lastRow, ...) {
    vcat(modelinfo$verbose, "...preparing permuted error matrices")
    
    if (firstRow > lastRow)
        stop("first row can't be greater than the last row")
    prows <- firstRow:lastRow
    nperms <- length(prows)
    if (nperms > ncol(p))
        stop("more rows requested than permutations")
    
    IHmat <- big.matrix(modelinfo$nobs^2, nperms+1, ...)
    IHmat[,1] <- as.vector(modelinfo$IH)
    
    l_ply(1:nperms, function(kk) {
        k <- prows[kk]
        prm <- p[,k]
        IHmat[,kk+1] <- as.vector(modelinfo$IH[prm,prm])
    }, .progress=modelinfo$progress)
    
    return(IHmat)
}

mdmr.prepare.pmat <- function(modelinfo, nvoxs, ...)
{
    vcat(modelinfo$verbose, "...preparing matrix of p-values")
    big.matrix(nvoxs, modelinfo$nfactors2perm, ...)
}

mdmr.prepare.fperms <- function(modelinfo, nperms, nvoxs, backingpath=NULL, ...)
{
    vcat(modelinfo$verbose, "...preparing matrix with all the Pseudo-F statistics")
    n <- modelinfo$nfactors2perm
    fpnames <- modelinfo$factor2perm.names
    lapply(1:n, function(i) {
        if (is.null(backingpath)) {
            bm <- big.matrix(nperms+1, nvoxs, ...)
        } else {
            bm <- big.matrix(nperms+1, nvoxs, 
                             backingpath=backingpath, 
                             backingfile=sprintf("fperms_%s.bin", fpnames[i]), 
                             descriptorfile=sprintf("fperms_%s.desc", fpnames[i]), 
                             ...)                       
        }
        return(bm)
    })
}

mdmr_worker <- function(firstVox, lastVox, Gmat, H2mats, IHmat, df.Res, df.Exp, Pmat, Fperms) {
    # n = number of subjects
    # where H2mat has rows of H's (n^2 rows) and cols of # of permutations
    # where Gmat has rows of dmats (n^2 rows)  and cols of # of voxels
    # result has rows of # of permutations and cols of # of voxels
    
    inds <- firstVox:lastVox
    nperms <- ncol(IHmat)
    nvoxs <- length(inds)
    nterms <- length(H2mats)
    
    # get part of bigmatrix with Gower's centered subject distances
    Gmat.chunk <- deepcopy(Gmat, cols=inds, shared=FALSE)
    
    # get part of Pmat
    Pmat.chunk <- sub.big.matrix(Pmat, firstRow=firstVox, lastRow=lastVox)
    
    # compute error term
    error.variance <- big.matrix(nperms, nvoxs, type="double", shared=FALSE)
    dgemm(C=error.variance, A=IHmat, B=Gmat.chunk, TRANSA='t')
    
    explained.variance <- big.matrix(nperms, nvoxs, type="double", shared=FALSE)
    for (i in 1:nterms) {
        # compute explained variance
        dgemm(C=explained.variance, A=H2mats[[i]], B=Gmat.chunk, TRANSA='t')
        
        # explained-variance / error
        Fstats <- sub.big.matrix(Fperms[[i]], firstCol=firstVox, lastCol=lastVox)
        do.operator(explained.variance, error.variance, "/", Fstats)
        
        # adjust for degrees of freedom
        dscal(Y=Fstats, ALPHA=df.Res/df.Exp[i])
        
        # get pvals (TODO: convert below line to C++ code)
        #Pmat.chunk[,i] <- apply(Fstats[,], 2, function(x) sum(x >= x[1])/nperms)
        .Call("ComputePvalsMain", Fstats@address, Pmat.chunk@address, as.double(i), PACKAGE="connectir")
    }
    
    # cleanup
    rm(explained.variance)
    rm(error.variance)
    gc(FALSE)
    
    return(NULL)
}

# assume each column of x has been gower centered
#mdmr <- function(G, formula, model, 
#                 nperms=4999, factors2perm=NULL, superblocksize=length(voxs), 
#                 voxs=1:ncol(G), blocksize=nperms, 
#                 contr.unordered="contr.sum", contr.ordered="contr.poly", 
#                 max.iter=10, strata=NULL, 
#                 verbose=1, parallel=FALSE, 
#                 G.path=NULL, fperms.path=NULL, 
#                 type="double", shared=parallel) 
#{
#    inform <- verbose==2
#    verbose <- as.logical(verbose)
#    progress <- ifelse(verbose, "text", "none")
#    vcat(verbose, "MDMR Time")
#    
#    if (!is.data.frame(model))
#        stop("input model must be a data frame")
#    if (!is.big.matrix(G) || !is.shared(G))
#        stop("input Gower's matrix must be type big matrix and shared")
#    if (!is.null(fperms.path) && !file.exists(fperms.path))
#        stop("output path for Pseudo-F permutations does not exist")
#    if (is.filebacked(G) && is.null(G.path))
#        stop("G.path must be given when G is filebacked big matrix")
#    if (!is.filebacked(G) && !is.null(G.path))
#        warning("G.path will be ignored", immediate.=TRUE)
#    if (!is.null(G.path) && !file.exists(G.path))
#        stop("output path for distance matrices does not exist")
#    
#    nvoxs <- length(voxs)
#    nsubs <- sqrt(nrow(G))  # don't want length(subs) here
#    if (nvoxs == 0) stop("# of voxels cannot be 0")
#    if (nperms == 0) stop("# of permutations cannot be 0")    
#    superblocks <- niftir.split.indices(1, nvoxs, by=superblocksize)
#    blocks <- niftir.split.indices(1, nperms, by=blocksize)
#    
#    # Prepare model matrix
#    modelinfo <- mdmr.prepare.model(formula, model, contr.unordered, contr.ordered, 
#                                    factors2perm, verbose)
#    vcat(verbose, "...will calculate p-values for the following factors: %s", 
#         paste(modelinfo$factor2perm.names, collapse=", "))
#    modelinfo$verbose <- FALSE
#    modelinfo$progress <- "none"
#    
#    # Permutation Business
#    p <- mdmr.prepare.permutations(modelinfo, nperms, strata, max.iter)
#    
#    # Create output matrices
#    Pmat <- mdmr.prepare.pmat(modelinfo, nvoxs, type=type, shared=TRUE)
#    if (!is.null(fperms.path)) {
#        save_fperms <- TRUE
#        Fperms <- mdmr.prepare.fperms(modelinfo, nperms, nvoxs, 
#                                      backingpath=fperms.path, 
#                                      type=type, shared=TRUE)
#    } else {
#        save_fperms <- FALSE
#        Fperms <- NULL
#    }
#    
#    dfRes <- as.double(modelinfo$df.Res)
#    dfExp <- as.double(modelinfo$df.Exp)
#    
#    vcat(verbose, "Computing MDMR across %i large blocks", 
#         superblocks$n)
#    vcat(verbose, "...with %i smaller blocks within each larger one",
#         blocks$n)
#    
#    # tmp_G <- deepcopy(G...sub_nvoxs)
#    # tmp_F <- nperms, sub_nvoxs
#    
#    for (si in 1:superblocks$n) {
#        vcat(verbose, "large block %i", si)
#        start.time <- Sys.time()
#        
#        firstVox <- superblocks$starts[si]
#        lastVox <- superblocks$ends[si]
#        sub_nvoxs <- lastVox - firstVox + 1
#        
#        tmp_Pmat <- mdmr.prepare.pmat(modelinfo, sub_nvoxs, init=1, 
#                                      type=type, shared=shared)
#        
#        vcat(verbose, "...copying distance matrices into memory")
#        tmpG <- deepcopy(x=G, cols=firstVox:lastVox, type=type, shared=shared)
#        if (is.filebacked(G))
#            G <- free.memory(G, G.path)
#        
#        vcat(verbose, "...preparing pseudo-F permutation matrices")
#        list_tmpFperms <- llply(1:blocks$n, function(bi) {
#            firstPerm <- as.double(blocks$starts[bi])
#            lastPerm <- as.double(blocks$ends[bi])
#            sub_nperms <- lastPerm - firstPerm + 1
#            tmpFperms <- mdmr.prepare.fperms(modelinfo, sub_nperms, sub_nvoxs, 
#                                             type=type, shared=shared)
#            return(tmpFperms)
#        }, .progress=progress, .inform=inform, .parallel=FALSE)
#        
#        afun <- function(bi) {
#            firstPerm <- as.double(blocks$starts[bi])
#            lastPerm <- as.double(blocks$ends[bi])
#            
#            tmpFperms <- list_tmpFperms[[bi]]
#            H2mats <- mdmr.prepare.permH2s(modelinfo, p, firstPerm, lastPerm, 
#                                           type=type, shared=FALSE)
#            IHmat <- mdmr.prepare.permIH(modelinfo, p, firstPerm, lastPerm, 
#                                         type=type, shared=FALSE)
#            
#            .Call("mdmr_worker", tmpG, tmpFperms, 
#                  tmp_Pmat, H2mats, IHmat, dfRes, dfExp, 
#                  PACKAGE="connectir")
#            
#            rm(H2mats, IHmat, tmpFperms)
#            invisible(gc(FALSE, TRUE))
#            
#            return(TRUE)
#        }
#        
#        vcat(verbose, "...almost MDMR time")
#        rets <- llply(1:blocks$n, afun, 
#                      .progress=progress, .inform=inform, .parallel=parallel)
#        
#        vcat(verbose, "...saving P-values")
#        sub_Pmat <- sub.big.matrix(Pmat, firstRow=firstVox, lastRow=lastVox)
#        deepcopy(x=tmp_Pmat, y=sub_Pmat)
#        rm(tmp_Pmat)
#        invisible(gc(FALSE, TRUE))
#        
#        if (save_fperms) {
#            n <- length(Fperms)
#            fpnames <- modelinfo$factor2perm.names
#            
#            Fperms <- llply(1:n, function(fi) {
#                vcat(verbose, "...saving permuted pseudo-F statistics for '%s'", 
#                     fpnames[fi])
#                aF <- Fperms[[fi]]
#                
#                # first row
#                tF <- list_tmpFperms[[1]][[fi]]
#                sF <- sub.big.matrix(aF, firstRow=1, lastRow=1, 
#                                    firstCol=firstVox, lastCol=lastVox, 
#                                    backingpath=fperms.path)
#                deepcopy(x=tF, rows=1, y=sF)
#                
#                # all other rows
#                l_ply(1:blocks$n, function(bi) {
#                    firstPerm <- blocks$starts[bi] + 1
#                    lastPerm <- blocks$ends[bi] + 1
#                    tF <- list_tmpFperms[[bi]][[fi]]
#                    sF <- sub.big.matrix(aF, firstRow=firstPerm, lastRow=lastPerm, 
#                                        firstCol=firstVox, lastCol=lastVox, 
#                                        backingpath=fperms.path)
#                    deepcopy(x=tF, rows=2:nrow(tF), y=sF)
#                }, .progress=progress)
#                
#                # clear memory a bit
#                flush(aF)
#                aF <- free.memory(aF, fperms.path)
#                
#                return(aF)
#            })
#        }
#        
#        rm(tmpG, list_tmpFperms)
#        invisible(gc(FALSE, TRUE))
#        
#        # how long?
#        end.time <- Sys.time()
#        time.total <- as.numeric(end.time-start.time, units="mins")
#        time.left <- time.total*(superblocks$n-si)
#        vcat(verbose, "...took %.1f minutes (%.1f minutes left)\n", 
#             time.total, time.left)
#    }
#    
#    # divide Pmat by # of perms + 1
#    .Call("mdmr_nmat_to_pmat", Pmat, as.double(nperms+1), PACKAGE="connectir")
#    
#    structure(
#        list(
#            nfactors=modelinfo$nfactors, 
#            modelinfo=modelinfo,
#            pvals=Pmat,
#            fstats=Fperms,
#            fpath=fperms.path, 
#            perms=p
#        ),
#        class="mdmr"
#    )
#}

# assume each column of x has been gower centered
# NOTE: here G must be file-backed 
mdmr.sge <- function(G.path, formula, model, 
                 nperms=4999, factors2perm=NULL, superblocksize=NULL, 
                 voxs=NULL, blocksize=nperms, 
                 contr.unordered="contr.sum", contr.ordered="contr.poly", 
                 max.iter=10, strata=NULL, 
                 verbose=1, parallel=FALSE, 
                 fperms.path=NULL, 
                 type="double", shared=parallel, 
                 forks=1, threads=1, njobs=NULL, ignore.proc.error=FALSE)
{
    library(Rsge)
    
    inform <- verbose==2
    verbose <- as.logical(verbose)
    progress <- ifelse(verbose, "text", "none")
    vcat(verbose, "MDMR Time")
    
    if (!is.data.frame(model))
        stop("input model must be a data frame")
    if (!is.null(fperms.path) && !file.exists(fperms.path))
        stop("output path for Pseudo-F permutations does not exist")
    if (!is.null(G.path) && !file.exists(G.path))
        stop("path for distance matrices must exist")
    
    G <- attach.big.matrix(G.path)
    G.dir <- dirname(G.path)
    if (is.null(voxs)) voxs <- 1:ncol(G)
    if (is.null(superblocksize)) superblocksize <- length(voxs)
    nvoxs <- length(voxs)
    nsubs <- sqrt(nrow(G))  # don't want length(subs) here
    superblocks <- niftir.split.indices(1, nvoxs, by=superblocksize)
    blocks <- niftir.split.indices(1, nperms, by=blocksize)
    sge.options(sge.user.options = sprintf("-S /bin/bash -pe mpi_smp %i", threads*forks))
    if (is.null(njobs)) {
        njobs <- 1
        use.sge <- FALSE
    } else {
        use.sge <- TRUE
    }
    
    # Prepare model matrix
    modelinfo <- mdmr.prepare.model(formula, model, contr.unordered, contr.ordered, 
                                    factors2perm, verbose)
    vcat(verbose, "...will calculate p-values for the following factors: %s", 
         paste(modelinfo$factor2perm.names, collapse=", "))
    modelinfo$verbose <- FALSE
    modelinfo$progress <- "none"
    
    # Permutation Business
    p <- mdmr.prepare.permutations(modelinfo, nperms, strata, max.iter)
    
    # Create output matrices
    if (!is.null(fperms.path)) {
        save_fperms <- TRUE
        Fperms <- mdmr.prepare.fperms(modelinfo, nperms, nvoxs, 
                                      backingpath=fperms.path, 
                                      type=type, shared=TRUE)
        fperms.desc <- lapply(Fperms, describe)
    } else {
        save_fperms <- FALSE
        Fperms <- NULL
        fperms.desc <- NULL
    }
    
    dfRes <- as.double(modelinfo$df.Res)
    dfExp <- as.double(modelinfo$df.Exp)
    
    vcat(verbose, "Computing MDMR across %i large blocks", 
         superblocks$n)
    vcat(verbose, "...with %i smaller blocks within each larger one",
         blocks$n)
    
    # tmp_G <- deepcopy(G...sub_nvoxs)
    # tmp_F <- nperms, sub_nvoxs
    
    
    
    # G but should just pass G.path and G.dir
    
    pvals.list <- sge.parLapply(1:superblocks$n, function(si) {
        vcat(verbose, "large block %i", si)
        start.time <- Sys.time()
        
        set_parallel_procs(forks, threads, verbose, ignore.proc.error)  
        
        firstVox <- superblocks$starts[si]
        lastVox <- superblocks$ends[si]
        sub_nvoxs <- lastVox - firstVox + 1
        
        tmp_Pmat <- mdmr.prepare.pmat(modelinfo, sub_nvoxs, init=1, 
                                      type=type, shared=shared)
        
        vcat(verbose, "...copying distance matrices into memory")
        G <- attach.big.matrix(G.path)
        tmpG <- deepcopy(x=G, cols=firstVox:lastVox, type=type, shared=shared)
        G <- free.memory(G, G.dir)
        
        vcat(verbose, "...preparing pseudo-F permutation matrices")
        list_tmpFperms <- llply(1:blocks$n, function(bi) {
            firstPerm <- as.double(blocks$starts[bi])
            lastPerm <- as.double(blocks$ends[bi])
            sub_nperms <- lastPerm - firstPerm + 1
            tmpFperms <- mdmr.prepare.fperms(modelinfo, sub_nperms, sub_nvoxs, 
                                             type=type, shared=shared)
            return(tmpFperms)
        }, .progress=progress, .inform=inform, .parallel=FALSE)
        
        afun <- function(bi) {
            firstPerm <- as.double(blocks$starts[bi])
            lastPerm <- as.double(blocks$ends[bi])
            
            tmpFperms <- list_tmpFperms[[bi]]
            H2mats <- mdmr.prepare.permH2s(modelinfo, p, firstPerm, lastPerm, 
                                           type=type, shared=FALSE)
            IHmat <- mdmr.prepare.permIH(modelinfo, p, firstPerm, lastPerm, 
                                         type=type, shared=FALSE)
            
            .Call("mdmr_worker", tmpG, tmpFperms, 
                  tmp_Pmat, H2mats, IHmat, dfRes, dfExp, 
                  PACKAGE="connectir")
            
            rm(H2mats, IHmat, tmpFperms)
            invisible(gc(FALSE, TRUE))
            
            return(TRUE)
        }
        
        vcat(verbose, "...almost MDMR time")
        rets <- llply(1:blocks$n, afun, 
                      .progress=progress, .inform=inform, .parallel=parallel)
        
        invisible(gc(FALSE, TRUE))
        
        if (save_fperms) {
            n <- length(fperms.desc)
            fpnames <- modelinfo$factor2perm.names
            
            l_ply(1:n, function(fi) {
                vcat(verbose, "...saving permuted pseudo-F statistics for '%s'", 
                     fpnames[fi])
                aF <- attach.big.matrix(fperms.desc[[fi]], backingpath=fperms.path)
                
                # first row
                tF <- list_tmpFperms[[1]][[fi]]
                sF <- sub.big.matrix(aF, firstRow=1, lastRow=1, 
                                    firstCol=firstVox, lastCol=lastVox, 
                                    backingpath=fperms.path)
                deepcopy(x=tF, rows=1, y=sF)
                
                # all other rows
                l_ply(1:blocks$n, function(bi) {
                    firstPerm <- blocks$starts[bi] + 1
                    lastPerm <- blocks$ends[bi] + 1
                    tF <- list_tmpFperms[[bi]][[fi]]
                    sF <- sub.big.matrix(aF, firstRow=firstPerm, lastRow=lastPerm, 
                                        firstCol=firstVox, lastCol=lastVox, 
                                        backingpath=fperms.path)
                    deepcopy(x=tF, rows=2:nrow(tF), y=sF)
                }, .progress=progress)
                
                # clear memory a bit
                flush(aF)
                aF <- free.memory(aF, fperms.path)
            })
        }
        
        rm(tmpG, list_tmpFperms)
        invisible(gc(FALSE, TRUE))
        
        # how long?
        end.time <- Sys.time()
        time.total <- as.numeric(end.time-start.time, units="mins")
        time.left <- time.total*(superblocks$n-si)
        vcat(verbose, "...took %.1f minutes (%.1f minutes left)\n", 
             time.total, time.left)
        
        # there is probably a smarter way to do this, then just copying it
        list(pmat=tmp_Pmat[,], firstVox=firstVox, lastVox=lastVox)
    }, debug=inform, trace=inform, 
    packages=c("connectir"), function.savelist=ls(), 
    cluster=use.sge, njobs=njobs)
    
    # combine pvalues for different voxels
    vcat(verbose, "...combining P-values")
    Pmat <- mdmr.prepare.pmat(modelinfo, nvoxs, type=type, shared=TRUE)
    for (l in pvals.list) {
        Pmat[l$firstVox:l$lastVox,] <- l$pmat
    }
    rm(pvals.list)
    
    # divide Pmat by # of perms + 1
    .Call("mdmr_nmat_to_pmat", Pmat, as.double(nperms+1), PACKAGE="connectir")
    
    # reload Fperms
    if (save_fperms) {
        rm(Fperms); invisible(gc(FALSE, TRUE))
        Fperms <- lapply(fperms.desc, function(desc) {
            attach.big.matrix(desc, backingpath=fperms.path)
        })
    }
    
    # free stuff
    rm(G); invisible(gc(FALSE, TRUE))
    
    # convert Pmat to matrix
    Pmat <- as.matrix(Pmat)
    if (is.vector(Pmat))
        Pmat <- matrix(Pmat, length(Pmat), 1)
    
    structure(
        list(
            nfactors=modelinfo$nfactors, 
            modelinfo=modelinfo, 
            pvals=as.matrix(Pmat),
            fstats=Fperms,
            fpath=fperms.path, 
            perms=p
        ),
        class="mdmr"
    )
}

#save_mdmr <- function(obj, voxs, sdir, mdir, formula, verbose=TRUE) {
#    vcat(verbose, "Saving MDMR Results")
#    
#    mdmr.output <- mdir
#    if (!file.exists(mdmr.output)) {
#        vcat(verbose, "...creating MDMR output directory")
#        dir.create(mdmr.output)
#    }
#    
#    mpath <- function(...) file.path(mdmr.output, ...)
#    
#    vcat(verbose, "...reading in brain mask")
#    seedfn <- file.path(sdir, "mask.nii.gz")
#    mask <- read.mask(seedfn)
#    total_nvoxs <- sum(mask)
#    header <- read.nifti.header(seedfn)
#    
#    # Model Stuff
#    modelinfo <- obj$modelinfo
#    modelinfo$perms <- obj$perms
#    ## save
#    vcat(verbose, "...saving model info")
#    save(modelinfo, file=mpath("modelinfo.rda"))
#    ## formula
#    cat(as.character(formula)[-c(1:2)], file=mpath("formula.txt"))
#    ## factors
#    vcat(verbose, "...saving factor names")
#    factornames <- modelinfo$factor2perm.names
#    nfactors <- length(factornames)
#    totext <- sprintf("# All Factors\n%s\n# Permuted Factors\n%s\n", 
#                paste(modelinfo$factor.names, collapse=", "), 
#                paste(factornames, collapse=", "))
#    cat(totext, file=mpath("factorinfo.txt"))
#    ## evs
#    vcat(verbose, "...saving evs")
#    evs <- modelinfo$rhs
#    write.csv(evs, file=mpath("evs.csv"))
#    ## models
#    # H2s <- modelinfo$H2s
#    ## residuals
#    vcat(verbose, "...saving residuals")
#    IH <- modelinfo$IH
#    write.table(IH[,], quote=F, row.names=F, col.names=F, file=mpath("residuals.2D"))
#    
#    vcat(verbose, "...saving p-values (1-p)")
#    Pmat <- obj$pvals
#    for (i in 1:nfactors) {
#        fn <- mpath(sprintf("pvals_%s.nii.gz", factornames[i]))
#        p <- vector("numeric", total_nvoxs)
#        p[voxs] <- 1 - Pmat[,i]
#        write.nifti(p, header, mask, outfile=fn, odt="float")
#    }
#    
#    vcat(verbose, "...saving FDR corrected p-values (1-p)")
#    CorrPmat <- big.matrix(nrow(Pmat), ncol(Pmat), type="double")
#    for (i in 1:nfactors) {
#        fn <- mpath(sprintf("fdr_pvals_%s.nii.gz", factornames[i]))
#        tmp <- p.adjust(Pmat[,i], "BH")
#        CorrPmat[,i] <- tmp
#        p <- vector("numeric", total_nvoxs)
#        p[voxs] <- 1 - tmp
#        write.nifti(p, header, mask, outfile=fn, odt="float")
#        rm(tmp)
#        gc(FALSE)
#    }
#    
#    vcat(verbose, "...saving z-statistics")
#    for (i in 1:nfactors) {
#        fn <- mpath(sprintf("zstats_%s.nii.gz", factornames[i]))
#        tmp <- qt(Pmat[,i], Inf, lower.tail=FALSE)
#        p <- vector("numeric", total_nvoxs)
#        p[voxs] <- tmp
#        write.nifti(p, header, mask, outfile=fn, odt="float")
#        rm(tmp)
#        gc(FALSE)
#    }
#    
#    vcat(verbose, "...saving FDR corrected z-statistics")
#    for (i in 1:nfactors) {
#        fn <- mpath(sprintf("fdr_zstats_%s.nii.gz", factornames[i]))
#        tmp <- qt(CorrPmat[,i], Inf, lower.tail=FALSE)
#        p <- vector("numeric", total_nvoxs)
#        p[voxs] <- tmp
#        write.nifti(p, header, mask, outfile=fn, odt="float")
#        rm(tmp)
#        gc(FALSE)
#    }
#    rm(CorrPmat)
#    gc(FALSE)
#    
#    if (!is.null(obj$clust)) {
#        vcat(verbose, "...saving cluster corrected p-values (1-p), etc")
#        for (i in 1:nfactors) {
#            fac <- factornames[i]
#            # cluster
#            fn <- mpath(sprintf("clust_%s.nii.gz", fac))
#            write.nifti(obj$clust[,i], header, mask, odt="int", outfile=fn)
#            # cluster p-values
#            fn <- mpath(sprintf("clust_pvals_%s.nii.gz", fac))
#            tmp <- (1-Pmat[,i]) * ((obj$clust[,i]>0)*1)
#            write.nifti(tmp, header, mask, odt="float", outfile=fn)
#        }
#    }
#}

