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


