test.computepvals <- function() {
    nr <- 5000; nc <- 500   
    x <- as.big.matrix(matrix(rnorm(nr*nc), nr, nc))
    
    system.time(y1 <- apply(x[,], 2, function(y) sum(y >= y[1])/nrow(x)))
    
    library(biganalytics)
    system.time(y2 <- apply(x, 2, function(y) sum(y >= y[1])/nrow(x)))
    
    system.time({
        y3 <- big.matrix(nc, 1, init=0, type="double")
        .Call("ComputePvalsMain", x@address, y3@address, as.double(1), PACKAGE="connectir")
    })
    
    all.equal(y1, y2)
    all.equal(y2, y3[,])
}

test.mdmr.prepare.model <- function() {
    model <- data.frame(group=factor(c(rep(c(1,2,3),30))), age=rnorm(90), random=rnorm(90))
    return(mdmr.prepare.model(.~group+age+random, model))
}

test.mdmr <- function() {
    nsubs <- 90
    nvoxs <- 120
    Dmat <- as.big.matrix(matrix(rnorm(nsubs^2 * nvoxs), nsubs^2, nvoxs))
    model <- data.frame(group=factor(c(rep(1:3,nsubs/3))), age=rnorm(nsubs), random=rnorm(nsubs))
    tmp1 <- mdmr(Dmat, .~group+age+random, model, nperms=999, factors.to.perm=c("group"), block.size=50)
    Xmodel <- model.matrix(~ group+age+random, data=model)
    pp <- test.philmdmr(Dmat, Xmodel, test.columns=2:3, nperm=999, permat=tmp1$perms[,])
    
    res <- all.equal(as.vector(tmp1$H2mats[[1]][,-1]), as.vector(pp$H2mat[,]))
    if (res) cat("H2mats are good\n")
    else cat("H2mats are NOT good\n")
    
    res <- all.equal(as.vector(tmp1$IHmat[,-1]), as.vector(pp$IHmat[,]))
    if (res) cat("IHmat are good\n")
    else cat("IHmat are NOT good\n")
    
    
}

test.philmdmr <- function(D.array, X, test.columns=ncol(X), permat=NULL, nperm=if (!is.null(permat)) ncol(permat) else 999, H2mat=NULL, IHmat=NULL, report.every=50, chunksize = 1000) {
    require(MASS)
  	n = nrow(X)
  	if (is.null(permat)) {
  		permat = matrix(NA, n, nperm)
  		for (k in 1:nperm) permat[ , k] = sample(n)
  	}
  	ndim = length(dim(D.array))
  	## warnings part
  	if (ndim!=2 & ndim!=3) stop("'D.array' must be either 3-way array of distance matrices, or a matrix whose columns are vectorized distance matrices")
    if (prod(dim(D.array)[-ndim]) != n^2) stop("Distance matrices in 'D.array' must be n x n, where n is the number of rows of 'X'")
    ## variable definition
    corrected.p = pvalue = F = rep(NA, dim(D.array)[ndim])
    maxpermF = NULL
    
    H = X %*% ginv(X)
    IH = diag(n) - H
    H2 = H - X[ , -test.columns] %*% ginv(X[ , -test.columns])
    if (is.null(H2mat)) {
        H2mat = IHmat = matrix(NA, n^2, nperm)
        for (k in 1:nperm) {
    	    if (k %% report.every==0) cat("Permutation", k, "\n")
    	    prm = permat[ , k]
     	    H2mat[ , k] = H2[prm, prm]
    	    IHmat[ , k] = IH[prm, prm]
      	}
    }
    
    ## divide into chunks
    dir.create("./temp")  # save chunk file into a temporary directory
    for (i in 1:ceiling(dim(D.array)[ndim]/chunksize))    {
    	cat("Chunk", i, "\n")
        ss = ((i-1)*chunksize+1): min(dim(D.array)[ndim], i*chunksize)
        if (ndim==2) D.mat = D.array[,ss]
        else if (ndim==3) D.mat = matrix(D.array[,,ss], ncol=length(ss))
        #Amat = -D.mat[,]*D.mat[,]/2
        Amat = Dmat[,]
        
        F[ss] = ((n-ncol(X)) / length(test.columns)) * as.vector(crossprod(Amat, as.vector(H2)) / crossprod(Amat, as.vector(IH)))
        permF.chunk = ((n-ncol(X)) / length(test.columns)) * crossprod(H2mat, Amat) / crossprod(IHmat, Amat)
        pvalue[ss]=apply(rbind(F[ss], permF.chunk), 2, function(v) sum(v[1]<=v)) / (nperm+1)
        maxpermF=apply(cbind(maxpermF,permF.chunk), 1, max)
        save(permF.chunk, file=paste("./temp/Chunk-", i, ".RData", sep=""))
    }

    for (i in 1:length(corrected.p))
        corrected.p[i] = (sum(F[i]<=maxpermF)+1)/(nperm+1)
    
    # recombining chunks of permF
    cat("Combining chunks of permuted-data statistics...\n")
    matrixsize.cutoff = 999*2000    # to decide whether to use a big.matrix obj to store the result
    if (nperm*dim(D.array)[ndim] <= matrixsize.cutoff)    permF = matrix(NA, nrow=nperm, ncol=dim(D.array)[ndim])
    else  {
        require(bigmemory)
        permF = filebacked.big.matrix(nperm, dim(D.array)[ndim], backingfile="permF.bin", descriptorfile="permF.desc")
    }
    for (i in 1:ceiling(dim(D.array)[ndim]/chunksize))    {
            load(paste("./temp/Chunk-", i, ".RData", sep=""))
            ss = ((i-1)*chunksize+1): min(dim(D.array)[ndim], i*chunksize)
            permF[,ss] = permF.chunk
        }
    return(list(F=F, permF=permF, maxpermF=maxpermF, pvalue=pvalue, corrected.p=corrected.p, fdr=p.adjust(pvalue, "BH"), H2mat=H2mat, IHmat=IHmat))
}
