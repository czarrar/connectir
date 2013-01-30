#new_dmat <- function(sdist, col, bp=NULL, asdist=FALSE) {
#    scol <- sub.big.matrix(sdist, firstCol=col, lastCol=col, backingpath=bp)
#    dcol <- matrix(0, nrow(scol), ncol(scol))
#    dcopy(X=scol, Y=dcol)
#    dim(dcol) <- rep(sqrt(nrow(scol)), 2)
#    if (asdist)
#        return(as.dist(dcol))
#    else
#        return(dcol)
#}

vox_base_cross <- function(FUN, funclist1, y, blocksize, 
                              funclist2=funclist1, 
                              verbose=TRUE, parallel=FALSE, shared=parallel, 
                              ztransform=FALSE, design_mat=NULL, 
                              recursive.unlist = TRUE, 
                              ...)
{
    vcat(verbose, "...checks")
    if (length(funclist1) != length(funclist2))
        stop("length mismatch between first and second set of functionals")
    if (length(funclist1) != length(y))
        stop("length mismatch between first set of functionals and labels")
    
    vcat(verbose, "...setup")
    to.regress <- ifelse(is.null(design_mat), FALSE, TRUE)
    nsubs <- length(funclist1)
    nvoxs1 <- ncol(funclist1[[1]])
    nvoxs2 <- ncol(funclist2[[1]])
    voxs1 <- 1:nvoxs1
    voxs2 <- 1:nvoxs2
    blocks <- niftir.split.indices(1, nvoxs1, by=blocksize)
    progress <- ifelse(verbose, "text", "none")
    
    gfun <- function(bi) {
        vcat(verbose, "...block %i", bi)
        
        first <- blocks$starts[bi]; last <- blocks$ends[bi]
        sub_voxs <- first:last
        sub_nvoxs <- length(sub_voxs)
                
        # correlate
        vcat(verbose, "....correlate")
        subs.cormaps <- vbca_batch3(funclist1, c(first, last), funclist2, 
                                    ztransform=ztransform, 
                                    type="double", shared=shared)
        
        # loop through each seed
        vcat(verbose, '....classify')
        sfun <- function(si) {
            # combine and transpose correlation maps for 1 seed
            seedCorMaps <- big.matrix(nsubs, nvoxs2, type="double", shared=FALSE)
            .Call("subdist_combine_and_trans_submaps", subs.cormaps, as.double(si), 
                  as.double(voxs2), seedCorMaps, PACKAGE="connectir")
            if (to.regress) {
                r_seedCorMaps <- big.matrix(nsubs, nvoxs2, type="double", shared=FALSE)
                qlm_residuals(seedCorMaps, design_mat, FALSE, r_seedCorMaps)
                rm(seedCorMaps)
                seedCorMaps <- r_seedCorMaps
            }
                        
            #X <- matrix(NA, nsubs, nvoxs2)
            #dcopy(X=seedCorMaps, Y=X)
            X <- seedCorMaps[,]
            
            # some classification?
            res <- FUN(X, y, ...)
            
            rm(seedCorMaps); gc(FALSE, TRUE)
            
            return(res)
        }
        
        res <- llply(1:sub_nvoxs, sfun, .parallel=parallel, .progress=progress, .inform=T)
                
        rm(subs.cormaps); gc(FALSE, TRUE)
        
        return(res)
    }
    
    vcat(verbose, "...%i blocks", blocks$n)
    res <- lapply(1:blocks$n, gfun)
    res <- unlist(res, recursive=recursive.unlist)
        
    return(res)
}

vox_svm_cross <- function(funclist1, y, blocksize, 
                          funclist2=funclist1, 
                          cross=10, kernel="linear", 
                          verbose=T, parallel=F, shared=parallel, 
                          ztransform=TRUE, 
                          ret.fit=FALSE, 
                          ...) 
{
    vcat(verbose, "SVM on Subject Distances")
    library(e1071)
        
    # TODO: test this function
    if (ret.fit) {
        FUN <- function(X, y, cross, kernel, ...) {
            fit <- svm(X, y, cross=cross, kernel=kernel, fitted=FALSE, ...)
            fit
        }
        ru <- FALSE
    } else {
        FUN <- function(X, y, cross, kernel, ...) {
            fit <- svm(X, y, cross=cross, kernel=kernel, fitted=FALSE, ...)
            res <- c(fit$tot.accuracy, fit$tot.MSE[[1]])
            res[!is.null(res)]
        }
        ru <- TRUE
    }
    
    tmp <- vox_base_cross(FUN, funclist1, y, blocksize, funclist2, verbose, 
                            parallel, shared, ztransform, recursive.unlist=ru,  
                            cross=cross, kernel=kernel, ...)
    
    return(tmp)
}

vox_glmnet_cross <- function(funclist1, y, blocksize, 
                             funclist2=funclist1, 
                             cross=10, family=NULL, standardize=T, alpha=0.5, 
                             verbose=T, parallel=F, shared=parallel, 
                             ztransform=TRUE, 
                             ret.fit=FALSE, 
                             ...)
{
    vcat(verbose, "GLMnet on Subject Distances")
    library(glmnet)
    
    if (is.null(family)) {
        if (is.factor(y)) {
            family <- ifelse(length(levels(y))==2, "binomial", "multinomial")
        } else {
            family <- "gaussian"
        }
    }
    
    if (family %in% c("binomial", "multinomial"))
        type.measure <- "class"
    else
        type.measure <- "deviance"
    
    # TODO: test this function
    if (ret.fit) {
        FUN <- function(X, y, cross, family, standardize, type.measure, alpha, ...) {
            fit <- cv.glmnet(X, y, family=family, standardize=standardize, 
                             nfolds=cross, type.measure=type.measure, alpha=alpha, ...)
            fit
        }
        ru <- FALSE
    } else {
        FUN <- function(X, y, cross, family, standardize, type.measure, alpha, ...) {
            fit <- cv.glmnet(X, y, family=family, standardize=standardize, 
                             nfolds=cross, type.measure=type.measure, alpha=alpha, ...)
            min(fit$cvm)
        }
        ru <- TRUE
    }
    
    tmp <- vox_base_cross(FUN, funclist1, y, blocksize, funclist2, verbose, 
                            parallel, shared, ztransform, recursive.unlist=ru,  
                            cross=cross, family=family, standardize=standardize, 
                            type.measure=type.measure, alpha=alpha, ...)
    
    return(tmp)
}
