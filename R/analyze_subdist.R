base_analyze_subdist <- function(FUN, Xs, y, memlimit=2, bpath=NULL, 
                                 verbose=T, parallel=F, 
                                 recursive.unlist = TRUE, 
                                 ...) 
{
    vcat(verbose, "...checks/setup")
    if (!is.big.matrix(Xs))
        stop("input Xs must be big matrix")
    if (!is.shared(Xs))
        stop("input Xs must be a SHARED big matrix")
    if (is.filebacked(Xs) && is.null(bpath))
        stop("must specify bpath when Xs input is file-backed")
    if (!is.null(bpath) && !file.exists(bpath))
        vstop("couldn't find backing path '%s'", bpath)
    nvoxs <- ncol(Xs)
    nsubs <- length(y)
    nr <- nrow(Xs)
    if (nsubs != sqrt(nr))
        stop("length mismatch between 'y' labels and rows of Xs")
    
    vcat(verbose, "...setting memory limit of %f GB", memlimit)
    blocksize <- get_sdist_analysis_limit(memlimit, Xs)
    vcat(verbose, "...setting block size to %i (out of %i voxels)", blocksize, nvoxs)
    
    blocks <- niftir.split.indices(1, nvoxs, by=blocksize)
    progress <- ifelse(verbose, "text", "none")
    
    if (is.filebacked(Xs)) {
        vcat(verbose, "...re-loading subject distances")
        bfile <- describe(Xs)@description$filename
        dfile <- sprintf("%s.desc", rmext(bfile))
        tmp_Xs <- attach.big.matrix(file.path(bpath, dfile))
    }
    
    dmat <- matrix(0, nsubs, nsubs)
    afun <- function(i) {
        X <- sub.big.matrix(Xs, firstCol=i, lastCol=i, backingpath=bpath)
        dcopy(N=nr, Y=dmat, X=X)
        FUN(dmat, y, ...)
    }
    
    vcat(verbose, "...testing 1 voxel")
    test <- afun(1)
    if (is.null(test) || is.na(test))
        stop("...test failed")
    
    vcat(verbose, "...%i blocks", blocks$n)
    results <- list()
    for (bi in 1:blocks$n) {
        vcat(verbose, "...block %i", bi)
        first <- blocks$starts[bi]; last <- blocks$ends[bi]
        res <- llply(first:last, afun, 
                     .progress=progress, .inform=verbose, .parallel=parallel)
        results[[bi]] <- res
        Xs <- free.memory(Xs, bpath)
        gc(FALSE, TRUE)
    }
    
    return(unlist(results, recursive=recursive.unlist))
}

glmnet_subdist_cross <- function(Xs, y, memlimit=2, bpath=NULL, 
                                 family=NULL, standardize=TRUE, cross=10, 
                                 verbose=T, parallel=F, ...) 
{
    vcat(verbose, "SVM on Subject Distances")
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
    FUN <- function(dmat, y, cross, family, standardize, type.measure, ...) {
        fit <- cv.glmnet(dmat, y, family=family, standardize=standardize, 
                         nfolds=cross, type.measure=type.measure, alpha=0.5, ...)
        #ret <- c(cvm.min=min(fit$cvm), cvm.mean=mean(fit$cvm), 
        #         nzero.min=fit$nzero[fit$lambda==fit$lambda.min], 
        #         nzero.mean=mean(fit$nzero))
        #ret <- switch(result, 
        #    1 = fit$cvm[fit$lambda==fit$lambda.min], 
        #    2 = mean(fit$cvm), 
        #    3 = fit$nzero[fit$lambda==fit$lambda.min], 
        #    4 = mean(fit$nzero),
        #    vstop("unrecognized result %i for glmnet", result)
        #)
        min(fit$cvm)
    }
    
    tmp <- base_analyze_subdist(FUN, Xs, y, memlimit, bpath, verbose, parallel, 
                                recursive.unlist = TRUE, 
                                cross=cross, standardize=standardize, family=family, 
                                type.measure=type.measure, ...)
    
    return(tmp)
}

svm_subdist_cross <- function(Xs, y, memlimit=2, bpath=NULL, 
                                cross=10, kernel="linear", 
                                verbose=T, parallel=F, ...) 
{
    vcat(verbose, "SVM on Subject Distances")
    library(e1071)
        
    # TODO: test this function
    FUN <- function(dmat, y, cross, kernel, ...) {
        fit <- svm(dmat, y, cross=cross, kernel=kernel, fitted=FALSE, ...)
        res <- c(fit$tot.accuracy, fit$tot.MSE[[1]])
        res[!is.null(res)]
    }
    
    tmp <- base_analyze_subdist(FUN, Xs, y, memlimit, bpath, verbose, parallel, 
                                cross=cross, kernel=kernel, ...)
    
    return(tmp)
}

kmeans_subdist_cross <- function(Xs, y, memlimit=2, bpath=NULL, 
                                    iter.max=200, nstart=20, 
                                    verbose=T, parallel=F, ...) 
{
    vcat(verbose, "K-means on Subject Distances")
    
    k <- length(unique(y))
    if (k != 2)
        stop("kmeans subdist can only be done right now with 2 groups")
    
    # TODO: test this function
    FUN <- function(dmat, y, k, iter.max, nstart, ...) {
        ks <- kmeans(dmat, k, iter.max, nstart, ...)$cluster
        acc <- sum(diag(table(ks, y)))/length(y)
        if (acc < 0.5)  acc <- 1 - acc
        return(acc)
    }
    
    tmp <- base_analyze_subdist(FUN, Xs, y, memlimit, bpath, verbose, parallel, 
                                k=k, iter.max=iter.max, nstart=nstart, ...)
    
    return(tmp)
}


# TODO: same thing as above but for kmeans
