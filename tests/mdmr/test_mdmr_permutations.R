library(connectir)
library(testthat)
library(stringr)

context("MDMR - Permutations")

create_init_data <- function(n, ngrps=2) {
    n <- 100
    formula <- ~ group*continuous
    model <- data.frame(group = factor(rep(c(1:ngrps), each=n/ngrps)), 
                        continuous = rnorm(n))
    factor.names <- c("group", "continuous", "group:continuous")
    
    return(list(
        formula = formula, 
        model = model, 
        factor.names = factor.names
    ))
}

create_rhs <- function(n=100, ngrps=2) {
    dat <- create_init_data(n, ngrps)
    rhs <- mdmr_model.rhs(dat$formula, dat$model)
    return(rhs)
}

test_that("permuted index works without a strata", {
    n <- 100

    set.seed(10)
    ref.perms <- sample(n)
    
    set.seed(10)
    comp.perms <- permuted.index(n, NULL)
    
    expect_that(ref.perms, equals(comp.perms))    
})

test_that("permuted index works with strata", {
    n <- 100
    strata <- rep(c("Group1", "Group2"), each=n/2)

    set.seed(10)
    inds <- 1:n
    ref.perms <- c(sample(inds[strata=='Group1'], n/2), 
                   sample(inds[strata=='Group2'], n/2))
    
    set.seed(10)
    comp.perms <- permuted.index(n, strata)
    
    expect_that(ref.perms, equals(comp.perms))
})

test_that("can detect significantly correlated permutation", {
    rhs <- create_rhs(100)
    qrhs <- mdmr_model.qr(rhs)
    rhs <- mdmr_model.rank(rhs, qrhs)
    grps <- attr(qrhs, "grps")
    
    pinds <- 1:100
    pinds[c(1,2)] <- c(2,1)
    r <- abs(mdmr_perms.cor_with_perm(rhs, grps, 2, pinds))
    expect_that(r, equals(1))
    
    pinds <- 1:100
    pinds[c(1,100)] <- c(100,1)
    r <- abs(mdmr_perms.cor_with_perm(rhs, grps, 2, pinds))
    expect_that(r>0.9, is_true())    
})

test_that("can gather permutations that aren't significantly correlated with the model", {
    rhs <- create_rhs(100)
    qrhs <- mdmr_model.qr(rhs)
    rhs <- mdmr_model.rank(rhs, qrhs)
    grps <- attr(qrhs, "grps")
    factors2perm <- attr(qrhs, "factors2perm")
    
    pmats <- mdmr_perms.gather_perms_for_factor(rhs, grps, factors2perm[1], 50)
    
    # Check correlation with original model
    H2orig <- mdmr_model.hat_matrix_2(rhs, grps, factors2perm[1])
    cors <- apply(pmats, 2, function(ps) {
        H2perm <- mdmr_model.hat_matrix_2(rhs, grps, factors2perm[1], ps)
        abs(cor(H2orig[,1], H2perm[,1]))
    })
    rthresh <- tanh(1.96/sqrt(nrow(rhs)-3))
    
    expect_that(all(cors<rthresh), is_true())
})

test_that("permutes multiple columns for one factor with H2", {
    rhs <- create_rhs(100, 4)
    qrhs <- mdmr_model.qr(rhs)
    rhs <- mdmr_model.rank(rhs, qrhs)
    grps <- attr(qrhs, "grps")
    factors2perm <- attr(qrhs, "factors2perm")
    ps <- sample(100)
    
    # Reference
    TOL <- 1e-07
    ## H
    Xj <- rhs
    Xj[,2:4] <- Xj[ps,2:4]
    qrX <- qr(Xj, tol=TOL)
    Q <- qr.Q(qrX)
    H <- tcrossprod(Q[,1:qrX$rank])
    ## H2
    Xj <- rhs[,c(1,5:8)]
    qrX <- qr(Xj, tol = TOL)
    Q <- qr.Q(qrX)
    ref.H2 <- H - tcrossprod(Q[, 1:qrX$rank])
    
    # Comparison
    comp.H2 <- mdmr_model.hat_matrix_2(rhs, grps, factors2perm[1], ps)
    
    expect_that(ref.H2, equals(comp.H2))
})

test_that("permutes multiple columns for one factor with IH", {
    rhs <- create_rhs(100, 4)
    qrhs <- mdmr_model.qr(rhs)
    rhs <- mdmr_model.rank(rhs, qrhs)
    grps <- attr(qrhs, "grps")
    factors2perm <- attr(qrhs, "factors2perm")
    ps <- sample(100)
    
    # Reference
    TOL <- 1e-07
    ## H
    Xj <- rhs
    Xj[,2:4] <- Xj[ps,2:4]
    qrX <- qr(Xj, tol=TOL)
    Q <- qr.Q(qrX)
    H <- tcrossprod(Q[,1:qrX$rank])
    ## H2
    ref.IH <- diag(nrow(rhs)) - H
    
    # Comparison
    comp.IH <- mdmr_model.hat_matrix_ih(rhs, grps, factors2perm[1], ps)
    
    expect_that(ref.IH, equals(comp.IH))
})

test_that("gather permuted H2s for a factor", {
    rhs <- create_rhs(100)
    qrhs <- mdmr_model.qr(rhs)
    rhs <- mdmr_model.rank(rhs, qrhs)
    grps <- attr(qrhs, "grps")
    factors2perm <- attr(qrhs, "factors2perm")
    
    pmat <- mdmr_perms.gather_perms_for_factor(rhs, grps, factors2perm[1], 100, 
                                               include.orig=TRUE)
    
    # Reference
    ref.H2perms <- apply(pmat, 2, function(ps) {
        H2perm <- mdmr_model.hat_matrix_2(rhs, grps, factors2perm[1], ps)
        as.vector(H2perm)
    })
    
    # Comparison
    comp.H2perms <- mdmr_perms.gather_H2perms_for_factor(rhs, grps, 
                                                         factors2perm[1], pmat)
    
    expect_that(ref.H2perms, is_equivalent_to(comp.H2perms[,]))
})

test_that("gather permuted IHs for a factor", {
    rhs <- create_rhs(100)
    qrhs <- mdmr_model.qr(rhs)
    rhs <- mdmr_model.rank(rhs, qrhs)
    grps <- attr(qrhs, "grps")
    factors2perm <- attr(qrhs, "factors2perm")
    
    pmat <- mdmr_perms.gather_perms_for_factor(rhs, grps, factors2perm[1], 100)
    
    # Reference
    ref.IHperms <- apply(pmat, 2, function(ps) {
        IHperm <- mdmr_model.hat_matrix_ih(rhs, grps, factors2perm[1], ps)
        as.vector(IHperm)
    })
    
    # Comparison
    comp.IHperms <- mdmr_perms.gather_IHperms_for_factor(rhs, grps, 
                                                         factors2perm[1], pmat)
    
    expect_that(ref.IHperms, is_equivalent_to(comp.IHperms[,]))
})

test_that("gather permuted H2s", {
    rhs <- create_rhs(100)
    qrhs <- mdmr_model.qr(rhs)
    rhs <- mdmr_model.rank(rhs, qrhs)
    grps <- attr(qrhs, "grps")
    factors2perm <- attr(qrhs, "factors2perm")
    nfactors <- length(factors2perm)
    
    list.perms <- mdmr_perms.gather_perms(rhs, qrhs, 100, include.orig=TRUE)
    
    # Reference
    ref.H2s <- lapply(1:nfactors, function(fi) {
        perms <- list.perms[[fi]]
        apply(perms, 2, function(ps) {
            H2perm <- mdmr_model.hat_matrix_2(rhs, grps, factors2perm[fi], ps)
            as.vector(H2perm)
        })
    })
    
    # Comparison
    comp.H2s <- mdmr_perms.gather_H2s(rhs, qrhs, list.perms)
    
    for (fi in 1:nfactors)
        expect_that(ref.H2s[[fi]], is_equivalent_to(comp.H2s[[fi]][,]))
})

test_that("gather permuted IHs", {
    rhs <- create_rhs(100)
    qrhs <- mdmr_model.qr(rhs)
    rhs <- mdmr_model.rank(rhs, qrhs)
    grps <- attr(qrhs, "grps")
    factors2perm <- attr(qrhs, "factors2perm")
    nfactors <- length(factors2perm)
    
    list.perms <- mdmr_perms.gather_perms(rhs, qrhs, 100, include.orig=TRUE)
    
    # Reference
    ref.IHs <- lapply(1:nfactors, function(fi) {
        perms <- list.perms[[fi]]
        apply(perms, 2, function(ps) {
            IHperm <- mdmr_model.hat_matrix_ih(rhs, grps, factors2perm[fi], ps)
            as.vector(IHperm)
        })
    })
        
    # Comparison
    comp.IHs <- mdmr_perms.gather_IHs(rhs, qrhs, list.perms)
    
    for (fi in 1:nfactors)
        expect_that(ref.IHs[[fi]], is_equivalent_to(comp.IHs[[fi]][,]))
})
