vbca <- function(bigmat, cols, ztransform=FALSE, outmat=NULL, ...) {
    A <- deepcopy(bigmat, cols=cols)
#    A <- sub.big.matrix(bigmat, firstCol=firstCol, lastCol=lastCol)
	B <- bigmat
	if (is.null(outmat))
	    C <- big.matrix(ncol(A), ncol(B), init=0, type="double", ...)
	else if (ncol(A)==nrow(outmat) && ncol(B)==ncol(outmat))
	    C <- outmat
	else
	    stop("dimensions of outmat are out of whack")
	ALPHA <- 1/(nrow(B) - 1)
	dgemm(C=C, A=A, B=B, TRANSA='t', ALPHA=ALPHA)
	if (ztransform)
	    atanh(C)
	invisible(C)
}

vbca_batch <- function(subs.bigmats, cols, ztransform=FALSE, ...) {
    nsubjects <- length(subs.bigmats)
    lapply(1:nsubjects, function(i) vbca(subs.bigmats[[i]], cols, ztransform, ...))
}

# Inverse Covariance / glasso functions

pcor <- function(x, y=NULL, d=1) {
	oc <- cov(x, y)
	ic <- -1*solve(oc)
	tmp <- sqrt(diag(ic)) %*% matrix(1,1,n)
	r <- ic / tmp / t(tmp)
	diag(r) <- d
	return(r)
}

icov <- function(x, y=NULL, lamda=100, d=1) {
    oc <- cov(x, y)
    r <- norm_glasso(oc, lamda, d)
    return(r)
}

# oc = covariance matrix (or time-series if cov=TRUE)
norm_glasso <- function(oc, lamda=100, d=1) {
    library(glasso)

	n <- ncol(oc)
    
    # sparse inverse covariance matrix
    ic <- -1*glasso(oc/mean(diag(oc)), lamda/1000)$wi
    
    # use diaganoal to get normalized coefficients
    tmp <- sqrt(abs(diag(ic))) %*% matrix(1,1,n)
    r <- ic / tmp / t(tmp)
    
    # set diagonal
    diag(r) <- d
    
    return(r) 
}

## NEW CODE USING ARMADILLO

big_cor <- function(x, y=NULL, z=NULL, byrows=FALSE, 
                    x_firstCol=1, x_lastCol=ncol(x), 
                    y_firstCol=NULL, y_lastCol=NULL, 
                    z_firstCol=1, z_lastCol=NULL, 
                    glasso = FALSE, lamda=100, 
                    ...) 
{
    # Setup 'y'
    if (is.null(y)) {
        y <- x
        if (is.null(y_firstCol))
            y_firstCol <- x_firstCol
        if (is.null(y_lastCol))
            y_lastCol <- x_lastCol
    }
    if (is.null(y_firstCol))
        y_firstCol <- 1
    if (is.null(y_lastCol))
        y_lastCol <- ncol(y)
    
    # Setup/Check 'z'
    if (is.null(z)) {   # Create output matrix
        if (byrows) {
            nrow <- nrow(x)
            ncol <- nrow(y)
        } else {
            nrow <- x_lastCol - x_firstCol + 1
            ncol <- y_lastCol - y_firstCol + 1
        }
        z <- big.matrix(nrow, ncol, ...)
        if (z_firstCol!=1 || !is.null(z_lastCol))
            warning("user specifications of z_firstCol and z_lastCol not used (must give z)", 
                    immediate.=TRUE)
        z_firstCol <- 1; z_lastCol <- ncol
    } else if (is.null(z_lastCol)) {
        z_lastCol <- ncol(z)
    }
    
    # Check all inputs
    if (!is.big.matrix(x) || !is.big.matrix(y) || !is.big.matrix(z))
        stop("inputs must be big matrices")
    
    # How to perform the crossproduct?
    if (byrows) {
        method <- "big_tcor"
    } else {
        method <- "big_cor"
    }
    
    # Correlate!
    .Call(method, x, y, z, 
            as.double(x_firstCol), as.double(x_lastCol), 
            as.double(y_firstCol), as.double(y_lastCol), 
            as.double(z_firstCol), as.double(z_lastCol), 
            PACKAGE = "connectir")
    
    # Inverse Covariance
    if (glasso)
        z[,] <- norm_glasso(z[,], lamda=lamda)
    
    invisible(z)
}

big_ztransform <- function(x) {
    if (!is.big.matrix(x))
        stop("input must be a big matrix")
    .Call("big_ztransform", x, PACKAGE="connectir")
}

# incols => c(firstCol, lastCol)
# outcols => c(firstCol, lastCol)
vbca2 <- function(inmat, incols=c(1,ncol(inmat)), ztransform=FALSE, 
                  outmat=NULL, outcols=c(1,incols[2]), ...) 
{
	if (length(incols)!=2 || length(outcols)!=2)
	    stop("incols and outcols must be a vector of length 2: c(firstCol, lastCol)")
	
	# TODO: figure out better way to handle this
	if (is.null(outmat))
	    z_lastCol = NULL
	else
	    z_lastCol = outcols[2]
	
	outmat <- big_cor(x=inmat, z=outmat, 
	                  x_firstCol=incols[1], x_lastCol=incols[2], 
	                  y_firstCol=1, y_lastCol=ncol(inmat), 
	                  z_firstCol=outcols[1], z_lastCol=z_lastCol, 
	                  ...)
	if (ztransform) atanh(outmat)
	
	invisible(outmat)
}

vbca_batch2 <- function(subs.bigmats, cols, ztransform=FALSE, ...) {
    nsubjects <- length(subs.bigmats)
    lapply(subs.bigmats, function(x) vbca2(x, cols, ztransform, ...))
}

# TODO:
# allow code to do correlation between two different time-series objects
# for now allow that will look at every single voxel
# but for the second input can be some ROI type thing...

vbca3 <- function(inmat1, incols1=c(1,ncol(inmat1)), 
                  inmat2=NULL, incols2=NULL, 
                  ztransform=FALSE, 
                  outmat=NULL, outcols=c(1,incols2[2]), 
                  ...)
{
    if (is.null(inmat2) && !is.null(incols2))
        stop("cannot specify incols2 without specifying inmat2")
    if (is.null(inmat2))
        inmat2 <- inmat1
    if (is.null(incols2))
        incols2 <- c(1,ncol(inmat2))
    
    if (length(incols1)!=2 || length(incols2)!=2 || length(outcols)!=2)
        stop("incols and outcols must be a vector of length 2: c(firstCol, lastCol)")
    
    if (is.null(outmat))
        z_lastCol <- NULL
    else
        z_lastCol <- outcols[2]
    
    outmat <- big_cor(x=inmat1, y=inmat2, z=outmat, 
	                  x_firstCol=incols1[1], x_lastCol=incols1[2], 
	                  y_firstCol=incols2[1], y_lastCol=incols2[2], 
	                  z_firstCol=outcols[1], z_lastCol=z_lastCol, 
	                  ...)
	if (ztransform) atanh(outmat)
	
	invisible(outmat)
}

vbca_batch3 <- function(subs.bigmats1, cols1, subs.bigmats2=NULL, cols2=NULL, 
                        ztransform=FALSE, ...) 
{
    if (is.null(subs.bigmats2)) {
        lapply(subs.bigmats1, function(x) vbca2(x, cols1, ztransform, ...))
    } else {
        mapply(function(x, y) vbca3(x, cols1, y, cols2, ztransform, ...), 
                x=subs.bigmats1, y=subs.bigmats2)
    }
}
