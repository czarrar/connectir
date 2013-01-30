library(connectir)
library(testthat)
library(stringr)

context("Miscellaneous Functions used for Reading Functional Data")

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


context("detect_ftypes")

test_that("detect_ftypes returns the proper file-type", {
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


context("mask.read")

test_that("mask.read can read a NIFTI file", {
    file <- get_masklist("nifti")[1]
    ref <- read.mask(file)
    comp <- mask.read(file, "nifti")
    expect_that(ref, equals(comp))
})

test_that("mask.read can read a TEXT (SPACE) files", {
    file <- get_masklist("txt")[1]
    ref <- as.vector(as.matrix(read.table(file)))>0
    comp <- mask.read(file, "space")
    expect_that(ref, equals(comp))
})

test_that("mask.read fails for OTHER files", {
    file <- get_masklist("txt")[1]
    expect_that(mask.read(file, "other"), throws_error())
})


context("mask.auto")

test_that("mask.auto will correctly mask an example dataset", {
    mask_file <- get_masklist("nifti")[1]
    ref <- mask.read(mask_file, "nifti")
    
    nifti_file <- get_filelist("nifti")[1]
    reader <- gen_big_reader("nifti")
    comp <- mask.auto(nifti_file, reader)
    expect_that(ref, equals(comp))
})


context("overlap_masks")

test_that("overlap_masks based on masks", {
    ### NIFTI
    vcat(T, "For nifti")
    
    files <- get_masklist("nifti")
    
    ref.masks <- laply(files, mask.read, "nifti", .progress="text")
    ref.overlap <- colMeans(ref.masks)
    ref.mask <- ref.overlap == 1
    ref <- list(mask=ref.mask, overlap=ref.overlap, sub.masks=ref.masks)
    
    comp <- overlap_masks(files)
    
    expect_that(ref, equals(comp))
    
    
    ### TEXT
    vcat(T, "For text")
    
    files <- get_masklist("txt")
    
    ref.masks <- laply(files, mask.read, "space", .progress="text")
    ref.overlap <- colMeans(ref.masks)
    ref.mask <- ref.overlap == 1
    ref <- list(mask=ref.mask, overlap=ref.overlap, sub.masks=ref.masks)
    
    comp <- overlap_masks(files)
    
    expect_that(ref, equals(comp))
})

test_that("overlap_masks based on automasking", {
    ### NIFTI
    vcat(T, "For nifti")
    
    files <- get_filelist("nifti")
    
    reader <- gen_big_reader("nifti", type="double", shared=FALSE)
    ref.masks <- laply(files, mask.auto, reader, .progress="text")
    ref.overlap <- colMeans(ref.masks)
    ref.mask <- ref.overlap == 1
    ref <- list(mask=ref.mask, overlap=ref.overlap, sub.masks=ref.masks)
    
    comp <- overlap_masks(files, TRUE)
    
    expect_that(ref, equals(comp))
    
    
    ### TEXT
    vcat(T, "For text")
    
    files <- get_filelist("txt")
    
    reader <- gen_big_reader("space", type="double", shared=FALSE)
    ref.masks <- laply(files, mask.auto, reader, .progress="text")
    ref.overlap <- colMeans(ref.masks)
    ref.mask <- ref.overlap == 1
    ref <- list(mask=ref.mask, overlap=ref.overlap, sub.masks=ref.masks)
    
    comp <- overlap_masks(files, TRUE)
    
    expect_that(ref, equals(comp))
})



