library(connectir)
library(testthat)
library(stringr)

context("Computing Subject Distances via SGE")

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

context("Setup")

test_that("", {
    # Nifti
    files <- get_filelist("nifti")
    ref <- "nifti"
    comp <- detect_ftypes(files)
    expect_that(ref, equals(comp))
    
    # Afni 1D
    files <- get_filelist("1d")
    ref <- "space"
    comp <- detect_ftypes(files)
    expect_that(ref, equals(comp))
    
    # Txt
    files <- get_filelist("txt")
    ref <- "space"
    comp <- detect_ftypes(files)
    expect_that(ref, equals(comp))    
})


