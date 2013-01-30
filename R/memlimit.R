n2gb <- function(x) x*8/1024^3
gb2n <- function(x) x/8*1024^3

n2mb <- function(x) x*8/1024^2
mb2n <- function(x) x/8*1024^2

# opts => list(blocksize=0, memlimit=6, verbose=TRUE)
# return opts
get_subdist_memlimit <- function(opts, nsubs, nvoxs, subs.ntpts, nvoxs2=NULL) {
    vcat(opts$verbose, "Determining memory demands")

    nforks <- getDoParWorkers()
    mem_limit <- as.numeric(opts$memlimit)
    
    # RAM for functionals
    mem_used4func <- sum(sapply(subs.ntpts, function(x) n2gb(x*nvoxs)))
    if (!is.null(nvoxs2)) {
        vcat(opts$verbose, "...%.2f GB used for 1st set of functional data", mem_used4func)
        mem_used4func2 <- sum(sapply(subs.ntpts, function(x) n2gb(x*nvoxs2)))
        vcat(opts$verbose, "...%.2f GB used for 2nd set of functional data", mem_used4func2)
        mem_used4func <- mem_used4func + mem_used4func2
    } 
    vcat(opts$verbose, "...%.2f GB used for all functional data", mem_used4func)
    
    # RAM for 2 distance matrices
    mem_by_dists <- n2gb(nsubs^2 * 2)
    vcat(opts$verbose, "...%.1f MB used for 2 distance matrices", 
         mem_by_dists*1024)
    
    # RAM for 1 connectivity map
    if (is.null(nvoxs2))
        nvoxs2 <- nvoxs
    n2onemap <- nvoxs2 * nsubs
    mem_by_seeds <- n2gb(n2onemap)
    vcat(opts$verbose, "...%.1f MB used for 2 connectivity maps across %i subjects", 
         2*mem_by_seeds*1024, nsubs)
    
    # Temporary seed maps (varies with nforks)
    mem_tmp <- mem_by_seeds*nforks
    
    # memory varies based on chunking of # of seed voxels and # of distance matrices
    mem_limit <- as.numeric(opts$memlimit)
    mem_fixed <- mem_used4func + mem_tmp
    f <- function(d, s) {
        mem_limit - mem_fixed - d*mem_by_dists - s*mem_by_seeds
    }
    # note: for d: s <= d <= nvoxs
    
    # Auto-Set Block Size
    if (opts$blocksize == 0 || opts$superblocksize == 0) {
        # minimum amount of RAM needed
        min_mem_needed <- mem_fixed + 2*mem_by_dists + 2*nforks*mem_by_seeds
        
        # limit in RAM use
        if (mem_limit < min_mem_needed) {
            vstop(paste("You require at least %.2f GB of memory but are limited", 
                        "to %.2f GB. Please reset the --memlimit option."), 
                  min_mem_needed, mem_limit)
        } else {
            vcat(opts$verbose, paste("...memory limit is %.2f GB and a minimum", 
                                     "of %.2f GB is needed"), 
                 mem_limit, min_mem_needed)
        }
        
        # can use max of both?
        m <- f(nvoxs, nvoxs)
        if (m > 0) {
            opts$superblocksize <- nvoxs
            opts$blocksize <- nvoxs
        } else {
            vcat(opts$verbose, paste("...if you wanted to hold everything in memory",
                                     "you would need at least %.2f GB of RAM"), 
                               mem_limit-m)
            vcat(opts$verbose, "...autosetting superblocksize and blocksize")
            
            # super blocks (# of voxels in distance matrices to chunk)
            tmp <- c(2*nforks, floor(seq(0,1,by=0.002)[-1]*nvoxs))
            s.choices <- c(tmp, tmp[-1][tmp[-1] > 2*nforks])
            ds <- sapply(s.choices, function(s) {
                tryCatch(uniroot(f, c(s,nvoxs), s=s)$root, error=function(ex) NA)
            })
            w <- (length(ds) + 1) - which.min(rev(ds))
            d <- floor(ds[w])
            if (length(d) == 0 || d == 0) {
                # 2nd attempt
                ds <- sapply(s.choices, function(s) f(d=nvoxs, s=s))
                if (any(ds>0)) {
                    d <- nvoxs
                } else {
                    stop("Sh*%, you don't have enough RAM")
                }
            }
            opts$superblocksize <- d
            
            # regular blocks (# of seed voxels to chunk)
            s <- tryCatch(floor(uniroot(f, c(2,nvoxs), d=d)$root), 
                    error=function(ex) nvoxs)
            opts$blocksize <- s
        }
    }
    
    vcat(opts$verbose, "...setting super-block size to %i (out of %i)", 
         opts$superblocksize, nvoxs)
    vcat(opts$verbose, "...setting block size to %i (out of %i)", 
         opts$blocksize, nvoxs)
     
    # adjust for # of forks
    if (nforks > 1) {
        opts$blocksize <- floor(opts$blocksize/nforks)
        vcat(opts$verbose, "...adjusting block size to %i based on %i forks", 
             opts$blocksize, nforks)
    }

    # checks
    if (opts$blocksize < 1)
        stop("block size must be greater than 1 and less then # of voxels")
    if (opts$superblocksize < 1)
        stop("super block size is less than 1 and less than # of voxels")
    
    
    # calculate amount of memory that will be used
    d <- opts$superblocksize
    s <- opts$blocksize
    m <- f(d, s*nforks); mem_used <- mem_limit - m
    if (mem_used > mem_limit) {
        vstop("You require %.2f GB of memory but have a limit of %.2f GB", 
              mem_used, mem_limit)
    } else {
        vcat(opts$verbose, "...%.2f GB of RAM will be used", mem_used)
    }
        
    return(opts)
}

get_mdmr_memlimit <- function(opts, nsubs, nvoxs, nperms, nfactors) {
    vcat(opts$verbose, "Determining MDMR memory demands")
    
    len_dmat <- nsubs^2
    nforks <- getDoParWorkers()
    
    tmp <- n2gb(nfactors*nsubs*nperms*2)
    mem_perms <- ifelse(nsubs < 2^31, tmp, tmp/2)
    mem_pvals <- n2gb(nvoxs*nfactors*2*2)
    vcat(opts$verbose, "...%.1f MB used for permutation indices", mem_perms*1024)
    vcat(opts$verbose, "...%.1f MB used for p-values", mem_pvals*1024)
    
    mem_gmats <- n2gb(2*len_dmat)               # *nvoxs
    mem_fperms <- n2gb(2*2*nfactors*(nperms+1)) # *nvoxs
    vcat(opts$verbose, "...minimum of %.1f MB used for distance matrices", 
         2*mem_gmats*1024)
    vcat(opts$verbose, "...minimum of %.1f MB used for permuted pseudo-F stats", 
         2*mem_fperms*1024)
    mem_by_voxs <- mem_gmats + mem_fperms
    
    mem_h2s <- n2gb(nfactors*len_dmat)          # *nperms
    mem_ih <- n2gb(nfactors*len_dmat)           # *nperms
    vcat(opts$verbose, "...minimum of %.1f MB used for permuted model matrices", 
         2*mem_h2s*1024)
    vcat(opts$verbose, "...minimum of %.1f MB used for permuted error matrices", 
         2*mem_ih*1024)
    mem_by_perms <- mem_h2s + mem_ih
    
    mem_tmp <- n2gb(4*1*1)                       # *nvoxs*nperms
    vcat(opts$verbose, "...minimum of %.1f MB used for temporary matrices", 
         2*mem_tmp*1024)
    
    mem_fixed <- mem_perms + mem_pvals + mem_h2s + mem_ih
    
    # memory varies based on chunking of # of voxels and # of perms
    mem_limit <- as.numeric(opts$memlimit)
    f <- function(v, p) {
        mem_limit - mem_fixed - v*mem_by_voxs - (p+1)*mem_by_perms - v*(p+1)*mem_tmp
    }
    
    if (opts$blocksize == 0 || opts$superblocksize == 0) {
        # limit good?
        min_mem_needed <- mem_fixed + 2*mem_by_voxs + 2*mem_by_perms + 2*mem_tmp
        if (mem_limit < min_mem_needed) {
            vstop(paste("You require at least %.2f GB of memory but are limited", 
                        "to %.2f GB. Please reset the --memlimit option."), 
                  min_mem_needed, mem_limit)
        } else {
            vcat(opts$verbose, paste("...memory limit is %.2f GB and a minimum", 
                                     "of %.3f GB is needed"), 
                 mem_limit, min_mem_needed)
        }
        
        get_superblocksize <- function(v) {
            if (!is.null(opts$jobs) && v > round(nvoxs/opts$jobs)) {
                v <- ceiling(nvoxs/opts$jobs)
            }
            return(v)
        }
        
        # can use max of both?
        m <- f(nvoxs, nperms)
        if (m > 0) {
            opts$superblocksize <- get_superblocksize(nvoxs)
            opts$blocksize <- nperms
        } else {
            vcat(opts$verbose, paste("...if you wanted to hold everything in memory",
                                     "you would need at least %.2f GB of RAM"), 
                               mem_limit-m)
            vcat(opts$verbose, "...autosetting superblocksize and blocksize")
            
            # super blocks (# of voxels to chunk)
            m <- f(nvoxs, floor(nperms/2))
            if (m > 0) {
                v <- nvoxs
            } else {
                p.choices <- c(2*nforks, floor(seq(0,1,by=0.01)[-1]*nperms))
                vs <- sapply(p.choices, function(p) {
                    tryCatch(uniroot(f, c(2,nvoxs), p=p)$root, error=function(ex) NA)
                })
                if (all(is.na(vs))) {
                    # Just try to use all voxels
                    v <- nvoxs
                } else {
                    vs <- vs[!is.na(vs)]            
                    w <- (length(vs) + 1) - which.min(rev(vs))
                    v <- floor(vs[w])
                }
            }
            if (length(v) == 0 || v == 0) {
                stop("Sh*%, you don't have enough RAM (2)")
            } else {
                opts$superblocksize <- get_superblocksize(v)
            }
            
            # regular blocks (# of permutations to chunk)
            p <- tryCatch(floor(uniroot(f, c(2,nperms), v=v)$root), 
                    error=function(ex) nperms)
            opts$blocksize <- p
        }
    }
    
    vcat(opts$verbose, "...setting super-block size to %i (out of %i voxels)", 
         opts$superblocksize, nvoxs)
    vcat(opts$verbose, "...setting block size to %i (out of %i permutations)", 
         opts$blocksize, nperms)
    
    # adjust for # of forks
    if (nforks > 1) {
        opts$blocksize <- floor(opts$blocksize/nforks)
        vcat(opts$verbose, "...adjusting block size to %i based on %i forks", 
             opts$blocksize, nforks)
    }
    
    # checks
    if (opts$blocksize < 1)
        stop("block size is less than 1")
    if (opts$superblocksize < 1)
        stop("super block size is less than 1")
    
    # calculate amount of memory that will be used
    v <- opts$superblocksize
    p <- opts$blocksize
    m <- f(v, p*nforks); mem_used <- mem_limit - m
    if (mem_used > mem_limit) {
        vstop("You require %.2f GB of memory but have a limit of %.2f GB", 
              mem_used, mem_limit)
    } else {
        vcat(opts$verbose, "...%.2f GB of RAM will be used", mem_used)
    }
    
    return(opts)
}

.get_kendall_limit <- function(blocksize, mem_limit, nvoxs, subs.ntpts, nvoxs2=NULL, 
                              verbose=TRUE) 
{
    nsubs <- length(subs.ntpts)
    nforks <- getDoParWorkers()
    
    mem_func <- sum(sapply(subs.ntpts, function(x) n2gb(x*nvoxs)))
    if (!is.null(nvoxs2)) {
        vcat(verbose, "...%.2f GB used for 1st set of functional data", mem_func)
        mem_func2 <- sum(sapply(subs.ntpts, function(x) n2gb(x*nvoxs2)))
        vcat(verbose, "...%.2f GB used for 2nd set of functional data", mem_func2)
        mem_func <- mem_func + mem_func2
    } 
    vcat(verbose, "...%.2f GB used for all functional data", mem_func)
    mem_kendall <- n2gb(nvoxs)
    mem_fixed <- mem_func + mem_kendall
    
    if (is.null(nvoxs2))
        nvoxs2 <- nvoxs
    
    # for at least 1 seed
    mem_seedmap <- n2gb(nvoxs2*1*nforks)
    mem_cormaps <- n2gb(nvoxs2*nsubs*1*nforks)
    mem_by_seeds <- mem_seedmap + mem_cormaps
    
    mem_limit <- as.numeric(mem_limit)
    f <- function(s) {
        mem_limit - mem_fixed - s*mem_by_seeds
    }
    
    if (blocksize == 0) {
        # limit good?
        min_mem_needed <- mem_fixed + 2*mem_by_seeds
        if (mem_limit < min_mem_needed) {
            vstop(paste("You require at least %.2f GB of memory but are limited", 
                        "to %.2f GB. Please reset the --memlimit option."), 
                  min_mem_needed, mem_limit)
        } else {
            vcat(verbose, paste("...memory limit is %.2f GB and a minimum", 
                                     "of %.3f GB is needed"), 
                 mem_limit, min_mem_needed)
        }
        
        m <- f(nvoxs)
        if (m > 0) {
            blocksize <- nvoxs
        } else {
            vcat(verbose, paste("...if you wanted to hold everything in memory",
                                     "you would need at least %.2f GB of RAM"), 
                               mem_limit-m)
            vcat(verbose, "...autosetting blocksize")
            
            ss <- tryCatch(floor(uniroot(f, c(2,nvoxs))$root), error=function(ex) NA)
            w <- (length(ss) + 1) - which.min(rev(ss))
            s <- floor(ss[w])
            if (length(s) == 0 || s == 0) {
                stop("Sh*%, you don't have enough RAM")
            } else {
                blocksize <- s
            }
        }
    }
    
    vcat(verbose, "...setting block size to %i (out of %i voxels)", 
         blocksize, nvoxs)
    
    # adjust for # of forks
    if (nforks > 1) {
        blocksize <- floor(blocksize/nforks)
        vcat(verbose, "...adjusting block size to %i based on %i forks", 
             blocksize, nforks)
    }
    
    # checks
    if (blocksize < 1)
        stop("block size is less than 1")
    
    # calculate amount of memory that will be used
    s <- blocksize*nforks
    m <- f(s); mem_used <- mem_limit - m
    if (mem_used > mem_limit) {
        vstop("You require %.2f GB of memory but have a limit of %.2f GB", 
              mem_used, mem_limit)
    } else {
        vcat(verbose, "...%.2f GB of RAM will be used", mem_used)
    }
    
    return(blocksize)
}

# loop through a set of voxels
# calculate connectivity maps
# get fit for given seed map
# get contrasts (results)
.get_glm_limit <- function(blocksize, mem_limit, nvoxs, subs.ntpts, 
                          nvoxs2, verbose, nevs, ncons) 
{
    nsubs <- length(subs.ntpts)
    nforks <- getDoParWorkers()
    
    # fixed
    ## functionals
    mem_func <- sum(sapply(subs.ntpts, function(x) n2gb(x*nvoxs)))
    if (!is.null(nvoxs2)) {
        vcat(verbose, "...%.2f GB used for 1st set of functional data", mem_func)
        mem_func2 <- sum(sapply(subs.ntpts, function(x) n2gb(x*nvoxs2)))
        vcat(verbose, "...%.2f GB used for 2nd set of functional data", mem_func2)
        mem_func <- mem_func + mem_func2
    } 
    vcat(verbose, "...%.2f GB used for all functional data", mem_func)
    ## everything else
    if (is.null(nvoxs2))
        nvoxs2 <- nvoxs
    mem_evs <- n2gb(nevs*nsubs)
    mem_cons <- n2gb(ncons*nevs)
    mem_seedmap <- n2gb(nsubs*nvoxs2*nforks)
    mem_dd <- n2gb(nsubs*nsubs)
    mem_coefs1 <- n2gb(nevs*nvoxs2*nforks)       # qlm_fit
    mem_residuals <- n2gb(nsubs*nvoxs2*nforks)   # qlm_fit
    mem_mse <- n2gb(nvoxs2*nforks)               # qlm_fit
    mem_coefs2 <- n2gb(ncons*nvoxs2*nforks)      # qlm_contrast
    mem_stderrs <- n2gb(ncons*nvoxs2*nforks)     # qlm_contrast
    mem_tvals1 <- n2gb(ncons*nvoxs2*nforks)      # qlm_contrast
    mem_fixed <- mem_func + mem_evs + mem_cons + mem_seedmap + mem_dd + 
                 mem_coefs1 + mem_residuals + mem_mse + 
                 mem_coefs2 + mem_stderrs + mem_tvals1
    
    # for at least 1 seed
    mem_cormaps <- n2gb(nvoxs2*nsubs)
    mem_tvals2 <- n2gb(ncons*nvoxs2*2)
    mem_by_seeds <- mem_cormaps + mem_tvals2
    
    mem_limit <- as.numeric(mem_limit)
    f <- function(s) {
        mem_limit - mem_fixed - s*mem_by_seeds
    }
    
    if (blocksize == 0) {
        # limit good?
        min_mem_needed <- mem_fixed + 2*mem_by_seeds
        if (mem_limit < min_mem_needed) {
            vstop(paste("You require at least %.2f GB of memory but are limited", 
                        "to %.2f GB. Please reset the --memlimit option."), 
                  min_mem_needed, mem_limit)
        } else {
            vcat(verbose, paste("...memory limit is %.2f GB and a minimum", 
                                     "of %.3f GB is needed"), 
                 mem_limit, min_mem_needed)
        }
        
        m <- f(nvoxs)
        if (m > 0) {
            blocksize <- nvoxs
        } else {
            vcat(verbose, paste("...if you wanted to hold everything in memory",
                                     "you would need at least %.2f GB of RAM"), 
                               mem_limit-m)
            vcat(verbose, "...autosetting blocksize")
            
            ss <- tryCatch(floor(uniroot(f, c(2,nvoxs))$root), error=function(ex) NA)
            w <- (length(ss) + 1) - which.min(rev(ss))
            s <- floor(ss[w])
            if (length(s) == 0 || s == 0) {
                stop("Sh*%, you don't have enough RAM")
            } else {
                blocksize <- s
            }
        }
    }
    
    vcat(verbose, "...setting block size to %i (out of %i voxels)", 
         blocksize, nvoxs)
    
    # checks
    if (blocksize < 1)
        stop("block size is less than 1")
    
    # calculate amount of memory that will be used
    s <- blocksize
    m <- f(s); mem_used <- mem_limit - m
    if (mem_used > mem_limit) {
        vstop("You require %.2f GB of memory but have a limit of %.2f GB", 
              mem_used, mem_limit)
    } else {
        vcat(verbose, "...%.2f GB of RAM will be used", mem_used)
    }
    
    return(blocksize)
}


get_gcor_limit <- function(blocksize, mem_limit, nvoxs, ntpts, verbose=TRUE)
{
    
    nforks <- getDoParWorkers()
    
    mem_func <- n2gb(nvoxs*ntpts)
    mem_res <- n2gb(nvoxs)
    mem_fixed <- mem_func + mem_res
    
    mem_cormaps <- n2gb(nvoxs*1)
    
    mem_limit <- as.numeric(mem_limit)
    f <- function(s) {
        mem_limit - mem_fixed - s*mem_cormaps
    }
    
    if (blocksize == 0) {
        # limit good?
        min_mem_needed <- mem_fixed + 2*mem_cormaps
        if (mem_limit < min_mem_needed) {
            vstop(paste("You require at least %.2f GB of memory but are limited", 
                        "to %.2f GB. Please reset the --memlimit option."), 
                  min_mem_needed, mem_limit)
        } else {
            vcat(verbose, paste("...memory limit is %.2f GB and a minimum", 
                                     "of %.3f GB is needed"), 
                 mem_limit, min_mem_needed)
        }
        
        m <- f(nvoxs)
        if (m > 0) {
            blocksize <- nvoxs
        } else {
            vcat(verbose, paste("...if you wanted to hold everything in memory",
                                     "you would need at least %.2f GB of RAM"), 
                               mem_limit-m)
            vcat(verbose, "...autosetting blocksize")
            
            ss <- tryCatch(floor(uniroot(f, c(2,nvoxs))$root), error=function(ex) NA)
            w <- (length(ss) + 1) - which.min(rev(ss))
            s <- floor(ss[w])
            if (length(s) == 0 || s == 0) {
                stop("Sh*%, you don't have enough RAM")
            } else {
                blocksize <- s
            }
        }
    }
    
    vcat(verbose, "...setting block size to %i (out of %i voxels)", 
         blocksize, nvoxs)
    
    # adjust for # of forks
    if (nforks > 1) {
        blocksize <- floor(blocksize/nforks)
        vcat(verbose, "...adjusting block size to %i based on %i forks", 
             blocksize, nforks)
    }
    
    # checks
    if (blocksize < 1)
        stop("block size is less than 1")
    
    # calculate amount of memory that will be used
    s <- blocksize*nforks
    m <- f(s); mem_used <- mem_limit - m
    if (mem_used > mem_limit) {
        vstop("You require %.2f GB of memory but have a limit of %.2f GB", 
              mem_used, mem_limit)
    } else {
        vcat(verbose, "...%.2f GB of RAM will be used", mem_used)
    }
    
    return(blocksize)
}

get_sdist_analysis_limit <- function(mem_limit, Xs) {
    nvoxs <- ncol(Xs)
    nr <- nrow(Xs)
    mem_dmat <- n2gb(nr)
    mem_result <- n2gb(nvoxs)
    nforks <- getDoParWorkers()
    f <- function(v) mem_limit - mem_result - nforks*mem_dmat - v*nforks*mem_dmat
    if (f(nvoxs) > 0) {
        blocksize <- nvoxs
    } else {
        if (!is.filebacked(Xs))
            stop("not enough memory...you might want to give file-backed subject distances as input")        
        blocksize <- tryCatch(floor(uniroot(f, c(2, nvoxs))$root), error=function(ex) NA)
        if (is.na(blocksize))
            stop("not enough memory")
    }
    
    return(blocksize)
}
    