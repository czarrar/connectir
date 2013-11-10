library(connectir)
library(testthat)
library(stringr)

context("MDMR - Model")

create_init_data <- function(n) {
    n <- 100
    formula <- ~ group*continuous
    model <- data.frame(group = rep(c(0,1), each=n/2), 
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

test_that("right-hand matrix is properly created", {
    n <- 100
    dat <- create_init_data(n)
    
    # Reference - RHS creation
    rhs.frame <- model.frame(dat$formula, dat$model, drop.unused.levels = TRUE)
    op.c <- options()$contrasts
    options(contrasts = c("contr.sum", "contr.poly"))
    ref.rhs <- model.matrix(dat$formula, rhs.frame)
    options(contrasts = op.c)
    attr(ref.rhs, "factor.names") <- dat$factor.names
    
    # Comparison - RHS
    comp.rhs <- mdmr_model.rhs(dat$formula, dat$model)
    
    expect_that(ref.rhs, equals(comp.rhs))    
})

test_that("confirm that an elment with NA will throw error", {
    dat <- create_init_data(n)
    dat$model[10,2] <- NA
    expect_that(mdmr_model.rhs(dat$formula, dat$model), 
                throws_error("One of your factors has an empty element"))
})

test_that("take QR decomposition of RHS and add relevant attributes", {
    rhs <- create_rhs(100)
    
    # Reference
    qrhs <- qr(rhs)
    u.grps <- grps <- 0:3
    attr(qrhs, "grps") <- grps
    attr(qrhs, "u.grps") <- u.grps
    ## add factor names
    factor.names <- attr(rhs, "factor.names")[u.grps]
    attr(qrhs, "factor.names") <- factor.names
    ## add factors to permute
    factors2perm <- 2:length(u.grps)
    names(factors2perm) <- factor.names
    attr(qrhs, "factors2perm") <- factors2perm
    ## add columns to permute
    cols2perm <- lapply(factors2perm, function(fi) which(grps==u.grps[fi]))
    attr(qrhs, "cols2perm") <- cols2perm
    ## save qrhs
    ref.qrhs <- qrhs
    
    # Comparison
    comp.qrhs <- mdmr_model.qr(rhs)
    
    # Will check that QR is done properly
    # and that the following attributes are added to qrhs:
    # grps, u.grps, factor.names, factors2perm, and cols2perm
    expect_that(ref.qrhs, equals(comp.qrhs))
})

test_that("mdmr_model.qr will return an error if RHS has NULL for 'assign'", {
    rhs <- create_rhs(100)
    fn <- attr(rhs, "factor.names")
    rhs <- rhs[,]
    
    expect_that(mdmr_model.qr(rhs), throws_error())
})

test_that("mdmr_model.qr will return an error if RHS has NULL for 'factor.names'", {
    rhs <- create_rhs(100)
    attr(rhs, "assign") <- NULL
    
    expect_that(mdmr_model.qr(rhs), throws_error())
})

test_that("adjust a rank-deficient matrix", {
    dat <- create_init_data(n)
    dat$model <- cbind(dat$model, bad=dat$model$continuous)
    dat$formula <- ~ group*continuous + bad
    rhs <- mdmr_model.rhs(dat$formula, dat$model)
    
    qrhs <- mdmr_model.qr(rhs)
    
    # Reference
    ref.rhs <- rhs[,c(1,2,3,5)]
    
    # Comparison
    comp.rhs <- mdmr_model.rank(rhs, qrhs, throw_error=FALSE)
    
    expect_that(ref.rhs, equals(comp.rhs))
})

test_that("confirm that error will be thrown if RHS is rank-deficient", {
    rhs <- create_rhs(100)
    rhs <- cbind(rhs, bad=rhs[,2])
    qrhs <- mdmr_model.qr(rhs)
    
    expect_that(mdmr_model.rank(rhs, qrhs), throws_error())
})

test_that("generate H2s properly using 'rhs' permute option", {
    rhs <- create_rhs(100)
    qrhs <- mdmr_model.qr(rhs)
    rhs <- mdmr_model.rank(rhs, qrhs)
    
    nobs <- nrow(rhs)
    perms <- sample(1:nobs)
    
    # Reference
    ## basics
    TOL <- 1e-07
    grps <- attr(qrhs, "grps")
    u.grps <- attr(qrhs, "u.grps")
    ## model with variable of interest
    ref.H2s <- lapply(2:length(u.grps), function(j) {
        Xj <- rhs
        Xj[,grps %in% u.grps[j]] <- Xj[perms, grps %in% u.grps[j]]
        qrX <- qr(Xj, tol=TOL)
        Q <- qr.Q(qrX)
        H <- tcrossprod(Q)
        
        Xj <- rhs[, grps %in% u.grps[-j]]
        qrX <- qr(Xj, tol = TOL)
        Q <- qr.Q(qrX)
        H - tcrossprod(Q[,1:qrX$rank])
    })
    
    # Comparison
    comp.H2s <- lapply(2:length(u.grps), function(j) mdmr_model.hat_matrix_2(rhs, grps, j, perms, "rhs"))
    
    expect_that(ref.H2s, equals(comp.H2s))
})

test_that("generate H2 properly using 'hat' permute option", {
    rhs <- create_rhs(100)
    qrhs <- mdmr_model.qr(rhs)
    rhs <- mdmr_model.rank(rhs, qrhs)
    
    nobs <- nrow(rhs)
    perms <- sample(1:nobs)
    
    # Reference
    ## basics
    TOL <- 1e-07
    grps <- attr(qrhs, "grps")
    u.grps <- attr(qrhs, "u.grps")
    ## model with variable of interest
    ref.H2s <- lapply(2:length(u.grps), function(j) {
        Xj <- rhs
        qrX <- qr(Xj, tol=TOL)
        Q <- qr.Q(qrX)
        H <- tcrossprod(Q)
        
        Xj <- rhs[, grps %in% u.grps[-j]]
        qrX <- qr(Xj, tol = TOL)
        Q <- qr.Q(qrX)
        H2 <- H - tcrossprod(Q[,1:qrX$rank])
        
        H2[perms,perms]
    })
    
    # Comparison
    comp.H2s <- lapply(2:length(u.grps), function(j) mdmr_model.hat_matrix_2(rhs, grps, j, perms, "hat"))
    
    expect_that(ref.H2s, equals(comp.H2s))
})

test_that("generate H2 properly using 'hat_with_covariates' permute option", {
    rhs <- create_rhs(100)
    qrhs <- mdmr_model.qr(rhs)
    rhs <- mdmr_model.rank(rhs, qrhs)
    
    nobs <- nrow(rhs)
    perms <- sample(1:nobs)
    
    # Reference
    ## basics
    TOL <- 1e-07
    grps <- attr(qrhs, "grps")
    u.grps <- attr(qrhs, "u.grps")
    ## model with variable of interest
    ref.H2s <- lapply(2:length(u.grps), function(j) {
        Xj <- rhs
        qrX <- qr(Xj, tol=TOL)
        Q <- qr.Q(qrX)
        H <- tcrossprod(Q[,1:qrX$rank])
        
        Xj <- rhs[, grps %in% u.grps[-j]]
        qrX <- qr(Xj, tol = TOL)
        Q <- qr.Q(qrX)
        H2 <- H - tcrossprod(Q[,1:qrX$rank])
        
        Xj <- rhs[,grps %in% u.grps[j]]
        qrX <- qr(Xj, tol = TOL)
        Q <- qr.Q(qrX)
        H1 <- H - tcrossprod(Q[,1:qrX$rank])
        IH1 <- diag(nobs) - H1
        
        IH1 %*% H2[perms,perms] %*% IH1
    })
    
    # Comparison
    comp.H2s <- lapply(2:length(u.grps), function(j) mdmr_model.hat_matrix_2(rhs, grps, j, perms, "hat_with_covariates"))    
    expect_that(ref.H2s, equals(comp.H2s))
})

test_that("generation of hat matrix (H2) throws an error if o.inds is not valid", {
    rhs <- create_rhs(100)
    qrhs <- mdmr_model.qr(rhs)
    rhs <- mdmr_model.rank(rhs, qrhs)
    perms <- sample(1:nobs)[1:10]
    grps <- attr(qrhs, "grps")
    u.grps <- attr(qrhs, "u.grps")
    
    expect_that(mdmr_model.hat_matrix_2(rhs, grps, j, perms), throws_error())
})

test_that("generation of hat matrix (H2) throws an error if permute option is not valid", {
    rhs <- create_rhs(100)
    qrhs <- mdmr_model.qr(rhs)
    rhs <- mdmr_model.rank(rhs, qrhs)
    perms <- sample(1:nobs)
    grps <- attr(qrhs, "grps")
    u.grps <- attr(qrhs, "u.grps")
    
    expect_that(mdmr_model.hat_matrix_2(rhs, grps, j, perms, "junk"), throws_error())
})

test_that("generate IHs properly for rhs option", {
    rhs <- create_rhs(100)
    qrhs <- mdmr_model.qr(rhs)
    rhs <- mdmr_model.rank(rhs, qrhs)
    
    nobs <- nrow(rhs)
    perms <- sample(1:nobs)
    
    # Reference
    ## basics
    TOL <- 1e-07
    grps <- attr(qrhs, "grps")
    u.grps <- attr(qrhs, "u.grps")
    ## model with variable of interest
    ref.IH <- lapply(2:length(u.grps), function(j) {
        Xj <- rhs
        Xj[,grps %in% u.grps[j]] <- Xj[perms, grps %in% u.grps[j]]
        qrX <- qr(Xj, tol=TOL)
        Q <- qr.Q(qrX)
        H <- tcrossprod(Q[,1:qrX$rank])
        diag(nobs) - H
    })
    
    # Comparison
    comp.IH <- lapply(2:length(u.grps), function(j) mdmr_model.hat_matrix_ih(rhs, grps, j, perms, "rhs"))
    
    expect_that(ref.IH, equals(comp.IH))
})


test_that("generate IHs properly for 'hat' permute option", {
    rhs <- create_rhs(100)
    qrhs <- mdmr_model.qr(rhs)
    rhs <- mdmr_model.rank(rhs, qrhs)
    
    nobs <- nrow(rhs)
    perms <- sample(1:nobs)
    
    # Reference
    ## basics
    TOL <- 1e-07
    grps <- attr(qrhs, "grps")
    u.grps <- attr(qrhs, "u.grps")
    ## model with variable of interest
    ref.IH <- lapply(2:length(u.grps), function(j) {
        Xj <- rhs
        qrX <- qr(Xj, tol=TOL)
        Q <- qr.Q(qrX)
        H <- tcrossprod(Q[,1:qrX$rank])
        IH <- diag(nobs) - H
        IH[perms,perms]
    })
    
    # Comparison
    comp.IH <- lapply(2:length(u.grps), function(j) mdmr_model.hat_matrix_ih(rhs, grps, j, perms, "hat"))
    
    expect_that(ref.IH, equals(comp.IH))
})

test_that("generate IHs properly for 'hat_with_covariates' permute option", {
    rhs <- create_rhs(100)
    qrhs <- mdmr_model.qr(rhs)
    rhs <- mdmr_model.rank(rhs, qrhs)
    
    nobs <- nrow(rhs)
    perms <- sample(1:nobs)
    
    # Reference
    ## basics
    TOL <- 1e-07
    grps <- attr(qrhs, "grps")
    u.grps <- attr(qrhs, "u.grps")
    ## model with variable of interest
    ref.IH <- lapply(2:length(u.grps), function(j) {
        Xj <- rhs
        qrX <- qr(Xj, tol=TOL)
        Q <- qr.Q(qrX)
        H <- tcrossprod(Q)
        IH <- diag(nobs) - H
        
        Xj <- rhs[,grps %in% u.grps[j]]
        qrX <- qr(Xj, tol = TOL)
        Q <- qr.Q(qrX)
        H1 <- H - tcrossprod(Q[,1:qrX$rank])
        IH1 <- diag(nobs) - H1
        
        IH1 %*% IH[perms,perms] %*% IH1
    })
    
    # Comparison
    comp.IH <- lapply(2:length(u.grps), function(j) mdmr_model.hat_matrix_ih(rhs, grps, j, perms, "hat_with_covariates"))
    
    expect_that(ref.IH, equals(comp.IH))
})


test_that("confirm that generate hat matrices properly", {
    rhs <- create_rhs(100)
    qrhs <- mdmr_model.qr(rhs)
    rhs <- mdmr_model.rank(rhs, qrhs)
    
    grps <- attr(qrhs, "grps")
    u.grps <- attr(qrhs, "u.grps")
    
    factors2perm <- 2:length(u.grps)
    nobs <- nrow(rhs)
    perms <- sample(1:nobs)
    
    # Reference
    ref.H2s <- lapply(factors2perm, function(j) mdmr_model.hat_matrix_2(rhs, grps, j, perms))
    ref.IHs <- lapply(factors2perm, function(j) mdmr_model.hat_matrix_ih(rhs, grps, j, perms))
    
    # Comparison
    comp <- mdmr_model.hat_matrices(rhs, qrhs, perms)
    comp.H2s <- comp$H2s
    comp.IHs <- comp$IHs
    
    expect_that(ref.H2s, is_equivalent_to(comp.H2s))
    expect_that(ref.IHs, is_equivalent_to(comp.IHs))
})

test_that("Degrees of freedom are calculated right", {
    rhs <- create_rhs(100)
    qrhs <- mdmr_model.qr(rhs)
    rhs <- mdmr_model.rank(rhs, qrhs)
    
    # Reference
    n <- nrow(rhs)
    grps <- attr(qrhs, "grps")
    u.grps <- attr(qrhs, "u.grps")
    factors2perm <- attr(qrhs, "factors2perm")
    df.Res <- n - qrhs$rank
    df.Exp <- sapply(u.grps[factors2perm], function(i) sum(grps == i))
    ref.dfs <- list(res=df.Res, exp=df.Exp)
    
    # Comparison
    comp.dfs <- mdmr_model.df(qrhs)
    
    expect_that(ref.dfs, equals(comp.dfs))
})

test_that("prepare MDMR model", {
    dat <- create_init_data(100)
    formula <- dat$formula
    model <- dat$model
    
    # Reference
    rhs <- mdmr_model.rhs(formula, model)
    qrhs <- mdmr_model.qr(rhs)
    rhs <- mdmr_model.rank(rhs, qrhs)
    hats <- mdmr_model.hat_matrices(rhs, qrhs)
    dfs <- mdmr_model.df(qrhs)
    ref.ret <- list(
        formula = formula, 
        model = model, 
        rhs = rhs, 
        qrhs = qrhs, 
        H2s = hats$H2s, 
        IHs = hats$IHs, 
        df.res = dfs$res, 
        df.exp = dfs$exp
    )
    
    # Comparison
    comp.ret <- mdmr_model(formula, model)
    
    expect_that(ref.ret, equals(comp.ret))
})


