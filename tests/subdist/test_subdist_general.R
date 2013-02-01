library(connectir)
library(testthat)
library(stringr)


context("General setup functions towards computing subject distances")

get_filelist <- function(type)
{
    basedir <- system.file("data", package="connectir")
    
    if (type == "nifti") {
        ext <- "nii.gz"
    } else if (type == "1d") {
        ext <- "1D"
    } else if (type == "txt") {
        ext <- "txt"
    } else {
        vstop("unrecognized type %s", type)
    }
    
    files <- file.path(basedir, sprintf("test_func%02i.%s", 1:9, ext))
    
    return(files)
}

get_masklist <- function(type) {
    basedir <- system.file("data", package="connectir")
    
    if (type == "nifti") {
        ext <- "nii.gz"
    } else if (type == "txt") {
        ext <- "txt"
    } else {
        vstop("unrecognized type %s", type)
    }
    
    files <- file.path(basedir, sprintf("test_mask%02i.%s", 1:9, ext))
    
    return(files)
}

get_inlist <- function(type) {
    nifti_files <- get_filelist(type)
    mask_files <- get_masklist(type)
    grp_mask_file <- "data/test_mask_grp.nii.gz"
        
    inlist <- load_funcs.prepare(nifti_files)
    inlist <- load_funcs.mask(inlist, automask=TRUE, 
                              subject.masks=mask_files, 
                              group.mask=grp_mask_file)
    
    return(inlist)
}

get_dists <- function(inlist1, inlist2) {
    
    dists_list <- create_subdist(outdir, infiles1, mask1, infiles2, mask2, opts, shared=parallel_forks)
}


context("...reading regressors")

test_that("read_regressors reads good for matrices", {
    # Create sample design matrix
    mat <- matrix(rnorm(10*3), 10, 3)
    colnames(mat) <- c("one", "two", "three")
    dfile <- tempfile(fileext=".txt")
    write.table(mat, file=dfile, row.names=F)
    
    # Reference
    ref <- mat
    
    # Comparison
    comp <- subdist.read_regressors(dfile)
    
    expect_that(ref, equals(comp))
    
    # Delete design matrix
    file.remove(dfile)
})

test_that("read_regressors reads good for data frames", {
    # Create sample design matrix
    df <- data.frame(group=rep(c("ADHD", "Control"), each=5), age=rnorm(10), iq=rnorm(10)*100)
    dfile <- tempfile(fileext=".txt")
    write.table(df, file=dfile, row.names=F)
    f <- ~ group + age + iq
    
    # Reference
    new_df <- read.table(dfile, header=T)
    rhs.frame <- model.frame(f, new_df, drop.unused.levels = TRUE)
    rhs <- model.matrix(f, rhs.frame)
    ref <- as.matrix(rhs[,])
    
    # Comparison
    comp <- subdist.read_regressors(dfile, f)
    
    expect_that(ref, equals(comp))
    
    # Delete design matrix
    file.remove(dfile)
})

test_that("read_regressors fails for rank deficient matrix", {
    # Create sample design matrix
    mat <- matrix(rnorm(10*3), 10, 3)
    mat <- cbind(mat, mat[,3])
    colnames(mat) <- c("one", "two", "three", "four")
    dfile <- tempfile(fileext=".txt")
    write.table(mat, file=dfile, row.names=F)
    
    # Test
    expect_that(subdist.read_regressors(dfile), throws_error())
    
    # Delete design matrix
    file.remove(dfile)    
})

test_that("read_regressors fails when any value is empty", {
    # Create sample design matrix
    mat <- matrix(rnorm(10*3), 10, 3)
    mat[4,2] <- NA
    colnames(mat) <- c("one", "two", "three")
    dfile <- tempfile(fileext=".txt")
    write.table(mat, file=dfile, row.names=F)
    
    # Test
    expect_that(subdist.read_regressors(dfile), throws_error())
})



context("...reading regressors")

test_that("memory_limit returns the size of inputs when supplied", {
    for (type in c("nifti", "txt")) {
        inlist <- get_inlist("nifti")
        
        # Reference
        blocksize <- 1
        superblocksize <- 2
        ref <- list(blocksize=blocksize, superblocksize=superblocksize)
        
        # Comparison
        comp <- subdist.memory_limit(4, blocksize, superblocksize, inlist, 
                                     nforks=1)
        
        expect_that(ref, expect_that(comp))
    }
})

test_that("memory_limit returns some size when inputs not supplied", {
    for (type in c("nifti", "txt")) {
        inlist <- get_inlist("nifti")
        
        test <- subdist.memory_limit(4, 0, 0, inlist)
        
        expect_that(test$blocksize==0, is_false())
        expect_that(test$superblocksize==0, is_false())
    }
})


