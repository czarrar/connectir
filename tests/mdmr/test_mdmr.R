library(connectir)
library(testthat)
library(stringr)

context("MDMR")

create_init_data <- function(n) {
    formula <- ~ group*continuous
    model <- data.frame(group = rep(c("A","B","C","D"), each=n/4), 
                        continuous = rnorm(n))
    factor.names <- c("group", "continuous", "group:continuous")
    
    return(list(
        formula = formula, 
        model = model, 
        factor.names = factor.names
    ))
}

create_rhs <- function(n=100) {
    dat <- create_init_data(n)
    rhs <- mdmr_model.rhs(dat$formula, dat$model)
    return(rhs)
}

create_Gs <- function(n, nvoxs) {
    I <- diag(n)
    ones <- matrix(1, nrow=n)
    C <- I - ones %*% t(ones)/n
    Gs.mat <- sapply(1:nvoxs, function(vi) {
        smat <- matrix(rnorm(nvoxs*n), nvoxs, n)
        dmat <- 1 - cor(smat)
        A <- -(dmat^2)/2
        G <- C %*% A %*% C
        as.vector(G)
    })
    as.big.matrix(Gs.mat)
}

test_that("MDMR worker computes the F-statistic appropriately", {
    nobs <- 100; nvoxs <- 10; nperms <- 100
    dat <- create_init_data(nobs)
    formula <- dat$formula
    model <- dat$model
    Gs <- create_Gs(nobs, nvoxs)
    
    modelinfo <- mdmr_model(formula, model)
    factor.names <- names(attr(modelinfo$qrhs, "factors2perm"))
    nfactors <- length(factor.names)
    
    list.perms <- mdmr_perms.gather_perms(modelinfo$rhs, modelinfo$qrhs, 
                                          nperms, include.orig=TRUE)
    
    # Comparison
    comp.Fperms <- lapply(1:nfactors, function(fi) big.matrix(nperms, nvoxs))
    mdmr.worker(modelinfo, Gs, list.perms, comp.Fperms)
    
    # Reference
    ref.Fperms <- lapply(1:nfactors, function(fi) big.matrix(nperms, nvoxs))
    list.H2s <- mdmr_perms.gather_H2s(modelinfo$rhs, modelinfo$qrhs, list.perms)
    list.IHs <- mdmr_perms.gather_IHs(modelinfo$rhs, modelinfo$qrhs, list.perms)
    for (fi in 1:nfactors) {
        MS.Exp <- crossprod(list.H2s[[fi]][,], Gs[,])/modelinfo$df.exp[fi]
        MS.Res <- crossprod(list.IHs[[fi]][,], Gs[,])/modelinfo$df.res
        ref.Fperms[[fi]][,] <- MS.Exp/MS.Res
    }
    
    for (fi in 1:nfactors)
        expect_that(ref.Fperms[[fi]][,], equals(comp.Fperms[[fi]][,]))
})

test_that("can convert p-values to f-statistics", {
    nfactors <- 3
    nperms <- 100
    nvoxs <- 100
    
    list.Fperms <- lapply(1:nfactors, function(fi) {
        mat <- matrix(rnorm(nperms*nvoxs), nperms, nvoxs)
        as.big.matrix(mat)
    })
    
    ref.pvals <- sapply(1:nfactors, function(fi) {
        fperms <- list.Fperms[[fi]]
        apply(fperms, 2, function(fs) {
            sum(fs>=fs[1])/nperms
        })
    })
    
    comp.pvals <- mdmr.fstats_to_pvals(list.Fperms)
    
    expect_that(ref.pvals, is_equivalent_to(comp.pvals))
})

test_that("MDMR works appropriately", {
    nobs <- 100; nvoxs <- 10; nperms <- 100
    dat <- create_init_data(nobs)
    formula <- dat$formula
    model <- dat$model
    Gs <- create_Gs(nobs, nvoxs)
    
    # Comparison
    comp.mdmr <- mdmr(Gs, formula, model, nperms=nperms, save.fperms=TRUE)
    comp.pvals <- comp.mdmr$pvals
    comp.Fperms <- comp.mdmr$fstats  
      
    # Reference
    modelinfo <- mdmr_model(formula, model)
    factor.names <- names(attr(modelinfo$qrhs, "factors2perm"))
    nfactors <- length(factor.names)
    ref.Fperms <- lapply(1:nfactors, function(fi) big.matrix(nperms+1, nvoxs))
    mdmr.worker(modelinfo, Gs, comp.mdmr$perms, ref.Fperms)
    ref.pvals <- mdmr.fstats_to_pvals(ref.Fperms)
    
    # Evaluate
    for (fi in 1:nfactors)
        expect_that(ref.Fperms[[fi]][,], equals(comp.Fperms[[fi]][,]))
    expect_that(ref.pvals[,], is_equivalent_to(comp.pvals[,]))
    expect_that(modelinfo, equals(comp.mdmr$modelinfo))
})
