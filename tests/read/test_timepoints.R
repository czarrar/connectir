library(connectir)
library(testthat)
library(stringr)

context("Batch Reading Number of Time-Points")

gen_inlist_nifti <- function() {
    # Path to files
    basedir <- system.file("data", package="connectir")
    files <- file.path(basedir, sprintf("test_func%02i.nii.gz", 1:9))
    
    # Type of files
    ftype <- detect_ftypes(files)
    
    # Function to read
    reader <- gen_big_reader(ftype, type="double", shared=TRUE)
    
    # Mask
    header <- read.nifti.header(files[1])
    mask <- rep(T, prod(header$dim[1:3]))
    
    # Output
    return(list(files=files, ftype=ftype, reader=reader, mask=mask))
}

gen_inlist_1D <- function() {
    # Path to files
    basedir <- system.file("data", package="connectir")
    files <- file.path(basedir, sprintf("test_func%02i.1D", 1:9))
    
    # Type of files
    ftype <- detect_ftypes(files)
    
    # Function to read
    reader <- gen_big_reader(ftype, type="double", shared=TRUE)
    
    # Mask
    tmp <- read.table(files[1])
    mask <- rep(T, ncol(tmp))
    
    # Output
    return(list(files=files, ftype=ftype, reader=reader, mask=mask))
}

gen_inlist_txt <- function() {
    # Path to files
    basedir <- system.file("data", package="connectir")
    files <- file.path(basedir, sprintf("test_func%02i.txt", 1:9))
    
    # Type of files
    ftype <- detect_ftypes(files)
    
    # Function to read
    reader <- gen_big_reader(ftype, type="double", shared=TRUE)
    
    # Mask
    tmp <- read.table(files[1])
    mask <- rep(T, ncol(tmp))
    
    # Output
    return(list(files=files, ftype=ftype, reader=reader, mask=mask))
}

test_that("Accurately determines # of time-points for 4D nifti files", {
    ref <- rep(5, 9)
    inlist <- gen_inlist_nifti()
    comp <- get_funclist_tpts(inlist)
    expect_that(ref, equals(ref))
})

test_that("Accurately determines # of time-points for 1D", {  
    ref <- rep(5, 9)
    inlist <- gen_inlist_1D()
    comp <- get_funclist_tpts(inlist)
    expect_that(ref, equals(ref))
})

test_that("Accurately determines # of time-points for text", {
    ref <- rep(5, 9)
    inlist <- gen_inlist_txt()
    comp <- get_funclist_tpts(inlist)
    expect_that(ref, equals(ref))
})
