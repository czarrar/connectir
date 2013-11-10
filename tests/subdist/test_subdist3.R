library(connectir)
library(testthat)
library(stringr)

context("Subject Distances v3")

get_inlist <- function(files, hdr=NULL, rois=NULL, ...) {    
    # Type of files
    ftype <- detect_ftypes(files)
    
    # Function to read
    reader <- gen_big_reader(ftype, type="double", ...)
    
    if (ftype == "nifti4d")
    
    # Mask
    header <- read.nifti.header(files[1])
    mask <- rep(T, prod(header$dim[1:3]))
    
    # Output
    return(list(files=files, ftype=ftype, reader=reader, mask=mask))
}

get_inlist_nifti <- function() {
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

get_inlist_1D <- function() {
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
    
    # ROIs
    rois <- 1:length(mask)
    
    # Header
    hdr <- read.nifti.header(file.path(basedir, "test_func01.nii.gz"))
    
    # Output
    return(list(
        files=files, 
        ftype=ftype, 
        reader=reader, 
        mask=mask, 
        rois=rois, 
        header=hdr
    ))
}

get_inlist_txt <- function() {
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
    
    # ROIs
    rois <- 1:length(mask)
    
    # Header
    hdr <- read.nifti.header(file.path(basedir, "test_func01.nii.gz"))
    
    # Output
    return(list(
        files=files, 
        ftype=ftype, 
        reader=reader, 
        mask=mask, 
        rois=rois, 
        header=hdr
    ))
}

load_inlist <- function(inlist, scale=FALSE, parallel=TRUE) {
    funclist <- load_and_mask_func_data2(inlist$files, inlist$reader, mask=inlist$mask, 
                                         verbose=TRUE, scale=scale, 
                                         type="double", shared=parallel)
    return(funclist)
}

set_memlimit <- function(inlist1, inlist2, memlimit=4) {
    # Ns for 1st list of inputs
    nsubs <- length(inlist1$files)
    nvoxs1 <- sum(inlist1$mask)
    ntpts1 <- get_funclist_tpts(inlist1)
    
    # Ns for 2nd list of inputs (note: assume same # of tpts)
    if (!is.null(inlist2)) {
        nvoxs2 <- sum(inlist2$mask)
        ntpts2 <- get_funclist_tpts(inlist2)
        for (i in 1:nsubs) {
            if (ntpts1[i] != ntpts2[i])
                vstop("subject #%i does not have the same # of timepoints for the 1st and 2nd functional datasets", i)
        }
    } else {
        nvoxs2 <- NULL
    }
    
    # List of options
    opts <- list(verbose=TRUE, memlimit=memlimit, blocksize=0, superblocksize=0)
    opts$"no-link-functionals" <- TRUE
    
    # New options with blocksize and superblocksize set
    opts <- get_subdist_memlimit(opts, nsubs, nvoxs1, ntpts1, nvoxs2)
    
    return(opts)
}

create_subdist_file <- function() {
    outdir <- tempdir()
    dists_list <- create_subdist(outdir, inlist1$files, inlist1$mask, inlist2$files, 
                                 inlist2$mask, opts, shared=T)
    return(dists_list)
}

# To test sge_wrapper3, will need to create a random set of data + save it


# In order to call the compute_subdist_sge_wrapper3, I need the following inputs:
## 1. inlist1 and inlist2 => get_inlist_nifti and get_inlist_1D
## 2. blocksize and superblocksize => set_memlimit
## 3. list.dists => create_subdist_file

# Step 1 (inlist)
inlist1 <- get_inlist_nifti()
inlist2 <- get_inlist_1D()

# Step 2 (sizes)
opts <- set_memlimit(inlist1, inlist2)

# Step 3 (dists)
list.dists <- create_subdist_file()


# 

# blocksize and superblocksize
# list.dists



# opts$ => verbose, memlimit, blocksize, superblocksize

opts <- get_subdist_memlimit(opts, nsubs, nvoxs1, ntpts1, nvoxs2)

