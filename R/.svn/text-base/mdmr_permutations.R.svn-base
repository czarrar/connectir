#' Generate permuted indices, H2s, and IHs for each factor
#' 
#' @author Zarrar Shehzad
#' @param rhs right-hand matrix
#' @param qrhs QR of right-hand matrix
#' @param nperms number of permutations
#' @param strata vector of length nobs that contains the indices that can be
#'               exchanged
#' @param max.iter maximum number of iterations for redoing a significantly 
#'                 correlated permutation
#' @param verbose boolean
#' @param list.perms use a previous list of permuted indices
#' @param include.orig boolean to include the original indices in returned 
#'        matrix. if true, this will count as the first permutation.
#' @param ... values given to big matrices that are created for H2 and IH
#' @return list of matrices including permutation indices, permuted H2s, 
#'         and permuted IHs
mdmr_perms <- function(rhs, qrhs, nperms, strata=NULL, max.iter=100, 
                        verbose=TRUE, list.perms=NULL, include.orig=FALSE, 
                        ...)
{
    if (is.null(list.perms)) {
        list.perms <- mdmr_perms.gather_perms(rhs, qrhs, nperms, strata, 
                                              max.iter, include.orig, verbose)
    }
    list.H2s <- mdmr_perms.gather_H2s(rhs, qrhs, list.perms, verbose, ...)
    list.IHs <- mdmr_perms.gather_IHs(rhs, qrhs, list.perms, verbose, ...)
        
    return(list(
        perms = list.perms, 
        H2s = list.H2s, 
        IHs = list.IHs
    ))
}


#' Generate permuted observation indices for each factor
#' 
#' @author Zarrar Shehzad
#' @param rhs right-hand matrix
#' @param qrhs QR of right-hand matrix
#' @param nperms number of permutations
#' @param strata vector of length nobs that contains the indices that can be
#'               exchanged
#' @param max.iter maximum number of iterations for redoing a significantly 
#'                 correlated permutation
#' @param include.orig boolean to include the original indices in returned 
#'        matrix. if true, this will count as the first permutation.
#' @param verbose boolean
#' @return list of matrices for each factor with permuted observation indices
mdmr_perms.gather_perms <- function(rhs, qrhs, nperms, strata=NULL, 
                                    max.iter=100, include.orig=FALSE, 
                                    verbose=TRUE)
{
    vcat(verbose, "Gathering permuted observation indices")
    
    grps <- attr(qrhs, "grps")
    factors2perm <- attr(qrhs, "factors2perm")
    factor.names <- names(factors2perm)
    
    list.perms <- lapply(1:length(factors2perm), function(i) {
        j <- factors2perm[i]
        name <- factor.names[i]
        vcat(verbose, "...for factor %s (#%i)", name, i)
        
        pmat <- mdmr_perms.gather_perms_for_factor(rhs, grps, j, nperms, 
                                                   strata, max.iter, 
                                                   include.orig, verbose)
        pmat <- mdmr_perms.to_integer(pmat)
        
        return(pmat)
    })
    names(list.perms) <- factor.names
    
    return(list.perms)
}


#' Generate permuted H2s for each factor
#' 
#' @author Zarrar Shehzad
#' @param rhs right-hand matrix
#' @param qrhs QR of right-hand matrix
#' @param list.perms List of matrices one for each factor. Each matrix is 
#'        nobs x nperms where each column contains the permuted indices.
#'        see \code{\link{mdmr_perms.gather_perms}}
#' @param verbose boolean
#' @param ... values given to big matrices that are created for H2
#' @return list of H2 matrices for each factor
mdmr_perms.gather_H2s <- function(rhs, qrhs, list.perms, verbose=TRUE, ...)
{
    vcat(verbose, "Generating permuted hat matrices (H2s)")
    
    grps <- attr(qrhs, "grps")
    factors2perm <- attr(qrhs, "factors2perm")
    factor.names <- names(factors2perm)
    if (any(names(list.perms) != factor.names))
        stop("names for list.perms doesn't match factor.names")
    
    list.H2s <- lapply(1:length(factors2perm), function(i) {        
        j <- factors2perm[i]
        name <- factor.names[i]
        
        vcat(verbose, "...for factor %s (#%i)", name, i)
        
        perms <- list.perms[[i]]
        H2mat <- mdmr_perms.gather_H2perms_for_factor(rhs, grps, j, perms, ...)
        
        return(H2mat)
    })
    names(list.H2s) <- factor.names
    
    return(list.H2s)
}


#' Generate permuted IHs for each factor
#' 
#' @author Zarrar Shehzad
#' @param rhs right-hand matrix
#' @param qrhs QR of right-hand matrix
#' @param list.perms List of matrices one for each factor. Each matrix is 
#'        nobs x nperms where each column contains the permuted indices.
#'        see \code{\link{mdmr_perms.gather_perms}}
#' @param verbose boolean
#' @param ... values given to big matrices that are created for H2
#' @return list of IH matrices for each factor
mdmr_perms.gather_IHs <- function(rhs, qrhs, list.perms, verbose=TRUE, ...)
{
    vcat(verbose, "Generating permuted hat matrices (IHs)")
    
    grps <- attr(qrhs, "grps")
    factors2perm <- attr(qrhs, "factors2perm")
    factor.names <- names(factors2perm)
    if (any(names(list.perms) != factor.names))
        stop("names for list.perms doesn't match factor.names")
    
    list.IHs <- lapply(1:length(factors2perm), function(i) {
        j <- factors2perm[i]
        name <- factor.names[i]
        vcat(verbose, "...for factor %s (#%i)", name, i)
        
        perms <- list.perms[[i]]
        IHmat <- mdmr_perms.gather_IHperms_for_factor(rhs, grps, j, perms, ...)
        
        return(IHmat)
    })
    names(list.IHs) <- factor.names
    
    return(list.IHs)
}


# BELOW: CODE TAKEN FROM vegan
#' Generates vector of permuted indices
#' @seealso vegan package with function of the same name
permuted.index <- function (n, strata) 
{
    if (missing(strata) || is.null(strata)) 
        out <- sample(n, n)
    else {
        if (length(strata) != n)
            stop("Length of strata must equal to n: ", strata)
        out <- 1:n
        inds <- names(table(strata))
        for (is in inds) {
            gr <- out[strata == is]
            if (length(gr) > 1) 
                out[gr] <- sample(gr, length(gr))
        }
    }
    out
}

# ABOVE: CODE DIRECTLY TAKEN FROM vegan

#' INTERNAL: Calculate correlation between original H2 and permuted H2
#' 
#' @author Zarrar Shehzad
#' @param rhs right-hand matrix
#' @param grps vector indicating the factor index for each regressor
#' @param f.ind factor index to examine
#' @param o.inds vector of indices for observations (e.g., for a permutation)
#' @param ... additional arguments for cor
#' @return correlation coefficient
mdmr_perms.cor_with_perm <- function(rhs, grps, f.ind, o.inds, ...)
{
    H2orig <- mdmr_model.hat_matrix_2(rhs, grps, f.ind)
    H2perm <- mdmr_model.hat_matrix_2(rhs, grps, f.ind, o.inds)
    cor(H2orig[,1], H2perm[,1], ...)
}

#' INTERNAL: Gather permutation indices for a given factor
#' 
#' @author Zarrar Shehzad
#' @param rhs right-hand matrix
#' @param grps vector indicating the factor index for each regressor
#' @param f.ind factor index to examine
#' @param nperms number of permutations
#' @param strata vector of length nobs that contains the indices that can be
#'               exchanged
#' @param max.iter maximum number of iterations for redoing a significantly 
#'                 correlated permutation
#' @param include.orig boolean to include the original indices in returned 
#'        matrix. if true, this will count as the first permutation.
#' @param verbose boolean
#' @return matrix of nobs x nperms
mdmr_perms.gather_perms_for_factor <- function(rhs, grps, f.ind, nperms, 
                                               strata=NULL, max.iter=100, 
                                               include.orig=FALSE, 
                                               verbose=TRUE)
{
    nobs <- nrow(rhs)
    # r value that is significant (p < 0.05)
    rthresh <- tanh(1.96/sqrt(nobs-3))
    progress <- ifelse(verbose, "text", "none")
    
    pmat <- laply(1:nperms, function(p) {
        for (i in 1:max.iter) {
            # permuted indices
            ps <- permuted.index(nobs, strata)
            # absolute correlation of original vs permuted model
            r <- abs(mdmr_perms.cor_with_perm(rhs, grps, f.ind, ps))
            # check significance of correlation
            if (r <= rthresh) {
                break
            } else if (i == max.iter) {
                msg <- sprintf(paste("Max iteration reached but permutation ", 
                                "still significantly correlated with model, ", 
                                "see column %i called %s", sep=""), j, 
                                names(factors2perm)[factors2perm == j])
                warning(msg)
            }
        }
        return(ps)
    }, .progress=progress)
    pmat <- t(pmat)
    
    if (include.orig) 
        pmat[,1] <- 1:nobs
    
    return(pmat)
}

mdmr_perms.to_integer <- function(pmat)
{
    nobs <- nrow(pmat)
    if (nobs < 2^31) {
        dims <- dim(pmat)
        pmat <- as.integer(pmat)
        dim(pmat) <- dims
    }
    return(pmat)
}

#' INTERNAL: Gather the permuted H2s for a given factor
#' 
#' @author Zarrar Shehzad
#' @param rhs right-hand matrix
#' @param grps vector indicating the factor index for each regressor
#' @param f.ind factor index to examine
#' @param pmat matrix of permutations indices of nobs x nperms
#' @param permute can be rhs, hat, or hat_with_covariates
#' @param ... values given to big.matrix to be returned
#' @return big.matrix of nobs^2 x nperms+1 with each column a permuted H2
mdmr_perms.gather_H2perms_for_factor <- function(rhs, grps, f.ind, pmat,  
                                                 permute="rhs", ...)
{
    nobs   <- nrow(pmat)
    nperms <- ncol(pmat)
    
    bigmat <- big.matrix(nobs^2, nperms, ...)
    
    for (p in 1:nperms) {
        H2 <- mdmr_model.hat_matrix_2(rhs, grps, f.ind, pmat[,p], permute=permute)
        bigmat[,p] <- as.vector(H2)
    }
    
    return(bigmat)
}

#' INTERNAL: Gather the permuted IHs for a given factor
#' 
#' @author Zarrar Shehzad
#' @param rhs right-hand matrix
#' @param grps vector indicating the factor index for each regressor
#' @param f.ind factor index to examine
#' @param pmat matrix of permutations indices of nobs x nperms
#' @param permute can be rhs, hat, or hat_with_covariates
#' @param ... values given to big.matrix to be returned
#' @return big.matrix of nobs^2 x nperms+1 with each column a permuted H2
mdmr_perms.gather_IHperms_for_factor <- function(rhs, grps, f.ind, pmat, 
                                                 permute="rhs", ...)
{
    nobs   <- nrow(pmat)
    nperms <- ncol(pmat)
    
    bigmat <- big.matrix(nobs^2, nperms, ...)
    
    for (p in 1:nperms) {
        IH <- mdmr_model.hat_matrix_ih(rhs, grps, f.ind, pmat[,p], permute)
        bigmat[,p] <- as.vector(IH)
    }
    
    return(bigmat)
}
