library(connectir)
library(testthat)
library(stringr)

context("qlm")

create_data <- function(n=100, k=2, m=3) {
    y <- matrix(rnorm(n*m), n, m)
    X <- cbind(rep(1, n), matrix(rnorm(n*k), n, k))
    return(list(y=y, X=X))
}

test_that("qlm dd works?", {
    dat <- create_data()
    
    ref <- diag(solve(t(dat$X) %*% dat$X))      # will return a vector
    comp <- qlm_dd(as.big.matrix(dat$X))                # will return 1 column
    
    expect_that(ref, is_equivalent_to(comp))
})

test_that("qlm fit works?", {
    dat <- create_data()
    
    ref <- lm(dat$y ~ dat$X - 1)
    comp <- qlm_fit(as.big.matrix(dat$y), as.big.matrix(dat$X))
    
    expect_that(as.vector(ref$coef), equals(as.vector(comp$coef[,])))
    
    expect_that(as.vector(ref$residuals), equals(as.vector(comp$residuals[,])))
    
    mse <- colSums(ref$residuals^2)/ref$df.residual
    expect_that(mse, is_equivalent_to(comp$mse[,]))
})

test_that("qlm rsquared works?", {
    # f <- ref.fit$fitted.values[,1]
    # mss <- sum((f - mean(f))^2)
    # r <- ref.fit$residuals[,1]
    # rss <- sum(r^2)
    # r2 <- mss/(mss+rss)
    # adj.r2 <- 1 - (1 - r2) * ((nrow(dat$X) - 1)/(nrow(dat$X) - ncol(dat$X)))
    
    dat <- create_data()
    
    ref.fit <- lm(dat$y ~ dat$X[,-1])
    comp.fit <- qlm_fit(as.big.matrix(dat$y), as.big.matrix(dat$X))
    
    tmp <- summary(ref.fit)
    ref1 <- as.vector(sapply(tmp, function(xx) xx$r.squared))
    ref2 <- as.vector(sapply(tmp, function(xx) xx$adj.r.squared))
    
    comp1 <- qlm_rsquared(as.big.matrix(dat$y), comp.fit, adjust=FALSE)
    comp2 <- qlm_rsquared(as.big.matrix(dat$y), comp.fit)
    
    expect_that(as.vector(ref1), equals(as.vector(comp1[,])))
    expect_that(as.vector(ref2), equals(as.vector(comp2[,])))
})


test_that("qlm residuals works?", {
    dat <- create_data(k=2)
    
    fit <- lm(dat$y ~ dat$X - 1)
    residuals <- qlm_residuals(as.big.matrix(dat$y), as.big.matrix(dat$X))
    
    ref <- as.matrix(fit$residuals)
    comp <- as.matrix(residuals)
    expect_that(ref, is_equivalent_to(comp))
})

test_that("qlm residuals works when adding back column means?", {
    dat <- create_data(k=2)
    
    fit <- lm(dat$y ~ dat$X - 1)
    residuals <- qlm_residuals(as.big.matrix(dat$y), as.big.matrix(dat$X), 
								add.mean=TRUE)
    
    ref <- as.matrix(fit$residuals)
	for (i in 1:ncol(ref)) ref[,i] <- ref[,i] + mean(dat$X[,i])
    comp <- as.matrix(residuals)
    expect_that(ref, is_equivalent_to(comp))
})

test_that("qlm contrasts works?", {
    dat <- create_data(k=2)
    
    fit <- qlm_fit(as.big.matrix(dat$y), as.big.matrix(dat$X))

    dd <- qlm_dd(as.big.matrix(dat$X))
    cons <- matrix(c(0,1,0,0,0,1), 2, 3, byrow=T)
    comp <- qlm_contrasts(fit, cons, dd)
    
    ref <- list()
    ref$coefficients <- cons %*% as.matrix(fit$coef)
    c_dd <- cons %*% dd
    ref$standard_errors <- sqrt(sapply(as.vector(fit$mse[,]), function(x) c_dd * x))
    ref$tvals <- ref$coefficients/ref$standard_errors
    
    expect_that(ref$coef, is_equivalent_to(as.matrix(comp$coef)))
    expect_that(ref$st, is_equivalent_to(as.matrix(comp$st)))
    expect_that(ref$tval, is_equivalent_to(as.matrix(comp$tval)))
})

test_that("sfun in vox_qlm works", {
    # 1. SETUP
    shared <- TRUE
    nvoxs <- 100; voxs <- 1:nvoxs
    sub_voxs <- 50:70; sub_nvoxs <- length(sub_voxs)
    nevs <- 2
    nsubs <- 10
    
    ys <- lapply(1:nsubs, function(x) matrix(rnorm(sub_nvoxs*nvoxs), sub_nvoxs, nvoxs))
    X <- cbind(rep(1,nsubs), matrix(rnorm(nsubs*nevs), nsubs, nevs))
    
    subs.cormaps <- lapply(ys, as.big.matrix)
    evs <- as.big.matrix(X)
    
    dd <- qlm_dd(evs)
    cons <- matrix(c(0,1,0,0,0,1), 2, 3, byrow=T)
    ncons <- nrow(cons)
    
    tmp_outmats1 <- lapply(1:ncons, function(x) {
                        big.matrix(nvoxs, sub_nvoxs, type="double", shared=shared)
                    })    
    tmp_outmats2 <- lapply(1:ncons, function(x) {
                        big.matrix(nvoxs, sub_nvoxs, type="double", shared=shared)
                    })    
    
    # 2. TEST
    sfun <- function(si) {
        # combine and transpose correlation maps for 1 seed
        seedCorMaps <- big.matrix(nsubs, nvoxs, type="double", shared=FALSE)
        .Call("subdist_combine_and_trans_submaps", subs.cormaps, as.double(si), 
              as.double(voxs), seedCorMaps, PACKAGE="connectir")
        
        # glm
        fit <- qlm_fit(seedCorMaps, evs, shared=FALSE)
        res <- qlm_contrasts(fit, cons, dd, shared=FALSE)
        
        for (ci in 1:ncons)
            tmp_outmats1[[ci]][,si] <- res$tvalues[ci,]
        
        rm(seedCorMaps, fit, res); gc(FALSE)

        return(NULL)
    }
    for (si in 1:sub_nvoxs)
        sfun(si)
    
    # 3. REFERENCE
    rfun <- function(si) {
        seedCorMaps <- t(sapply(subs.cormaps, function(x) x[si,]))
        fit <- lm(seedCorMaps ~ X)
        res <- summary(fit)
        for (ci in 1:ncons) {
            for (vi in 1:nvoxs) {
                tmp_outmats2[[ci]][vi,si] <- res[[vi]]$coefficients[ci+1,3]
            }
        }
        
        return(NULL)
    }
    for (si in 1:sub_nvoxs)
        rfun(si)
    
    for (ci in 1:ncons)
        expect_that(tmp_outmats1[[ci]][,], equals(tmp_outmats2[[ci]][,]))
})


test_that("gfun in vox_qlm works", {
    # 1. SETUP
    bp <- NULL; verbose <- TRUE
    ztransform <- T; shared <- TRUE; parallel <- FALSE; progress <- "text"
    ntpts <- 50
    nvoxs <- 100; voxs <- 1:nvoxs
    nevs <- 2
    nsubs <- 10
    blocksize <- 20
    blocks <- niftir.split.indices(1, nvoxs, by=blocksize)
    
    ys <- lapply(1:nsubs, function(x) scale(matrix(rnorm(ntpts*nvoxs), ntpts, nvoxs)))
    X <- cbind(rep(1,nsubs), matrix(rnorm(nsubs*nevs), nsubs, nevs))
    
    funclist <- lapply(ys, as.big.matrix)
    evs <- as.big.matrix(X)
    
    dd <- qlm_dd(evs)
    cons <- matrix(c(0,1,0,0,0,1), 2, 3, byrow=T)
    ncons <- nrow(cons)
    
    outmats1 <- lapply(1:ncons, function(ci) 
                        big.matrix(nvoxs, nvoxs, init=0, type="double", shared=shared))
    outmats2 <- lapply(1:ncons, function(ci) 
                        big.matrix(nvoxs, nvoxs, init=0, type="double", shared=shared))
    
    gfun <- function(bi) {
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
            sub_outmat <- sub.big.matrix(outmats1[[ci]], firstCol=first, lastCol=last)
            deepcopy(x=tmp_outmats[[ci]], y=sub_outmat)
            if (!is.null(bp)) {
                flush(sub_outmat); flush(outmats1[[ci]])
                outmats1[[ci]] <- free.memory(outmats1[[ci]], bp)
            }
        }
        
        rm(subs.cormaps, tmp_outmats); gc(FALSE)
    }
    
    vcat(verbose, "...%i blocks", blocks$n)
    for (bi in 1:blocks$n)
        gfun(bi)    
    
    rfun <- function(bi) {
        vcat(verbose, "...block %i", bi)
        
        first <- blocks$starts[bi]; last <- blocks$ends[bi]
        sub_voxs <- first:last
        sub_nvoxs <- length(sub_voxs)
        
        tmp_outmats <- lapply(1:ncons, function(x) {
                            big.matrix(nvoxs, sub_nvoxs, init=0, type="double", shared=shared)
                        })
        
        # correlate
        vcat(verbose, "....correlate")
        #subs.cormaps1 <- lapply(funclist, function(x) {
        #    cmap <- cor(x[,first:last], x[,])
        #    if (ztransform)
        #        cmap <- atanh(cmap)
        #    return(cmap)
        #})
        subs.cormaps <- vbca_batch2(funclist, c(first, last), 
                                    ztransform=ztransform, 
                                    type="double", shared=shared)
        
        # loop through each seed
        vcat(verbose, '....glm')
        sfun <- function(si) {
            seedCorMaps <- t(sapply(subs.cormaps, function(x) x[si,-voxs[sub_voxs[si]]]))
            fit <- lm(seedCorMaps ~ X)
            res <- summary(fit)
            for (ci in 1:ncons) {
                rvoxs <- voxs[-sub_voxs[si]]
                for (vi in 1:(nvoxs-1)) {
                    tmp_outmats[[ci]][rvoxs[vi],si] <- res[[vi]]$coefficients[ci+1,3]
                }
            }
            
            return(NULL)
        }
        
        llply(1:sub_nvoxs, sfun, .parallel=parallel, .progress=progress, .inform=T)
        
        # copy over temp data
        vcat(verbose, '....save')
        for (ci in 1:ncons) {
            sub_outmat <- sub.big.matrix(outmats2[[ci]], firstCol=first, lastCol=last)
            deepcopy(x=tmp_outmats[[ci]], y=sub_outmat)
            if (!is.null(bp)) {
                flush(sub_outmat); flush(outmats2[[ci]])
                outmats2[[ci]] <- free.memory(outmats2[[ci]], bp)
            }
        }
        
        rm(subs.cormaps, tmp_outmats); gc(FALSE)
    }
    
    vcat(verbose, "...%i blocks", blocks$n)
    for (bi in 1:blocks$n)
        rfun(bi)
    
    for (ci in 1:ncons)
        expect_that(outmats1[[ci]][,], equals(outmats2[[ci]][,]))
})


test_that("gfun in vox_qlm3 works", {
    # 1. SETUP
    bp <- NULL; verbose <- TRUE
    ztransform <- T; shared <- TRUE; parallel <- FALSE; progress <- "text"
    ntpts <- 50
    nvoxs1 <- 10; voxs1 <- 1:nvoxs1
    nvoxs2 <- 100; voxs2 <- 1:nvoxs2
    nevs <- 2
    nsubs <- 10
    blocksize <- 20
    blocks <- niftir.split.indices(1, nvoxs1, by=blocksize)
    
    ys1 <- lapply(1:nsubs, function(x) scale(matrix(rnorm(ntpts*nvoxs1), ntpts, nvoxs1)))
    ys2 <- lapply(1:nsubs, function(x) scale(matrix(rnorm(ntpts*nvoxs2), ntpts, nvoxs2)))
    X <- cbind(rep(1,nsubs), matrix(rnorm(nsubs*nevs), nsubs, nevs))
    
    funclist <- lapply(ys, as.big.matrix)
    evs <- as.big.matrix(X)
    
    dd <- qlm_dd(evs)
    cons <- matrix(c(0,1,0,0,0,1), 2, 3, byrow=T)
    ncons <- nrow(cons)
    
    outmats1 <- lapply(1:ncons, function(ci) 
                        big.matrix(nvoxs, nvoxs, init=0, type="double", shared=shared))
    outmats2 <- lapply(1:ncons, function(ci) 
                        big.matrix(nvoxs, nvoxs, init=0, type="double", shared=shared))
    
    gfun <- function(bi) {
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
            sub_outmat <- sub.big.matrix(outmats1[[ci]], firstCol=first, lastCol=last)
            deepcopy(x=tmp_outmats[[ci]], y=sub_outmat)
            if (!is.null(bp)) {
                flush(sub_outmat); flush(outmats1[[ci]])
                outmats1[[ci]] <- free.memory(outmats1[[ci]], bp)
            }
        }
        
        rm(subs.cormaps, tmp_outmats); gc(FALSE)
    }
    
    vcat(verbose, "...%i blocks", blocks$n)
    for (bi in 1:blocks$n)
        gfun(bi)    
    
    rfun <- function(bi) {
        vcat(verbose, "...block %i", bi)
        
        first <- blocks$starts[bi]; last <- blocks$ends[bi]
        sub_voxs <- first:last
        sub_nvoxs <- length(sub_voxs)
        
        tmp_outmats <- lapply(1:ncons, function(x) {
                            big.matrix(nvoxs, sub_nvoxs, init=0, type="double", shared=shared)
                        })
        
        # correlate
        vcat(verbose, "....correlate")
        #subs.cormaps1 <- lapply(funclist, function(x) {
        #    cmap <- cor(x[,first:last], x[,])
        #    if (ztransform)
        #        cmap <- atanh(cmap)
        #    return(cmap)
        #})
        subs.cormaps <- vbca_batch2(funclist, c(first, last), 
                                    ztransform=ztransform, 
                                    type="double", shared=shared)
        
        # loop through each seed
        vcat(verbose, '....glm')
        sfun <- function(si) {
            seedCorMaps <- t(sapply(subs.cormaps, function(x) x[si,-voxs[sub_voxs[si]]]))
            fit <- lm(seedCorMaps ~ X)
            res <- summary(fit)
            for (ci in 1:ncons) {
                rvoxs <- voxs[-sub_voxs[si]]
                for (vi in 1:(nvoxs-1)) {
                    tmp_outmats[[ci]][rvoxs[vi],si] <- res[[vi]]$coefficients[ci+1,3]
                }
            }
            
            return(NULL)
        }
        
        llply(1:sub_nvoxs, sfun, .parallel=parallel, .progress=progress, .inform=T)
        
        # copy over temp data
        vcat(verbose, '....save')
        for (ci in 1:ncons) {
            sub_outmat <- sub.big.matrix(outmats2[[ci]], firstCol=first, lastCol=last)
            deepcopy(x=tmp_outmats[[ci]], y=sub_outmat)
            if (!is.null(bp)) {
                flush(sub_outmat); flush(outmats2[[ci]])
                outmats2[[ci]] <- free.memory(outmats2[[ci]], bp)
            }
        }
        
        rm(subs.cormaps, tmp_outmats); gc(FALSE)
    }
    
    vcat(verbose, "...%i blocks", blocks$n)
    for (bi in 1:blocks$n)
        rfun(bi)
    
    for (ci in 1:ncons)
        expect_that(outmats1[[ci]][,], equals(outmats2[[ci]][,]))
})
