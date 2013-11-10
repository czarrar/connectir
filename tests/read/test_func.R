library(connectir)
library(testthat)
library(stringr)

context("Batch Reading Number of Functional Data")

gen_inlist_nifti <- function() 
{
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


context("Prepare: load_funcs.prepare")

test_that("load_funcs.prepare returns a list with files, ftype, and reader", {
    types <- c("nifti", "1d", "txt")
    for (type in types) {
        vcat(T, "type %s", type)
        files <- get_filelist(type)
    
        # Reference
        ftype <- detect_ftypes(files)
        reader <- gen_big_reader(ftype, type="double")
        ref <- list(files=files, ftype=ftype, reader=reader)
        class(ref) <- "inlist_prepare"
    
        # Comparison
        comp <- load_funcs.prepare(files)
    
        expect_that(ref, equals(comp))
    }
})




context("Mask: load_funcs.mask")

test_that("automask", {
    for (type in c("nifti", "txt")) {
        vcat(T, "For %s", type)
        
        files <- get_filelist(type)
    
        inlist <- load_funcs.prepare(files)
        auto <- overlap_masks(inlist$files, TRUE, TRUE)
        
        comp <- load_funcs.mask(inlist, automask=TRUE, detailed=TRUE)
        
        expect_that(auto, equals(comp$mask_details$auto))
    }
})

test_that("subject.masks compares with automask results", {
    for (type in c("nifti", "txt")) {
        vcat(T, "For %s", type)
        
        nifti_files <- get_filelist(type)
        mask_files <- get_masklist(type)
        
        # Reference
        inlist <- load_funcs.prepare(nifti_files)
        ref <- load_funcs.mask(inlist, automask=TRUE)
        
        # Comparison
        comp <- load_funcs.mask(inlist, automask=FALSE, subject.masks=mask_files)
        
        expect_that(ref, equals(comp))
    }
})

test_that("group", {
    for (type in c("nifti", "txt")) {
        vcat(T, "For %s", type)
        
        nifti_files <- get_filelist(type)
        mask_file <- "data/test_mask_grp.nii.gz"
        
        # Reference
        inlist <- load_funcs.prepare(nifti_files)
        ftype <- detect_ftypes(mask_file)
        ref <- mask.read(mask_file, ftype)
        
        # Comparison
        comp <- load_funcs.mask(inlist, group.mask=mask_file)
        
        expect_that(ref, equals(comp$mask))
    }
})

test_that("combine works with right data", {
    for (type in c("nifti", "txt")) {
        vcat(T, "For %s", type)
        
        nifti_files <- get_filelist(type)
        mask_files <- get_masklist(type)
        grp_mask_file <- "data/test_mask_grp.nii.gz" 
        
        # Reference
        inlist <- load_funcs.prepare(nifti_files)
        ref.1 <- load_funcs.mask(inlist, automask=TRUE)
        ref.2 <- load_funcs.mask(inlist, group.mask=grp_mask_file)
        ref <- ref.1$mask & ref.2$mask
        
        # Comparison
        comp <- load_funcs.mask(inlist, automask=TRUE, subject.masks=mask_files, 
                                group.mask=grp_mask_file)
        
        expect_that(ref, equals(comp$mask))
    }    
})

test_that("combine does not work with wrong data", {
    for (type in c("nifti", "txt")) {
        vcat(T, "For %s", type)
        
        nifti_files <- get_filelist(type)
        mask_files <- get_masklist(type)
        grp_mask_file <- "data/test_mask_grp.nii.gz"
        
        # Reference
        inlist <- load_funcs.prepare(nifti_files)
        ref.1 <- load_funcs.mask(inlist, automask=TRUE)
        ref.2 <- load_funcs.mask(inlist, subject.masks=mask_files)
        ref <- ref.1$mask & ref.2$mask
        
        # Comparison
        comp <- load_funcs.mask(inlist, automask=TRUE, subject.masks=mask_files, 
                                group.mask=grp_mask_file)
        
        expect_that(all(ref==comp$mask), is_false())
    }
})



context("load_funcs.read")

test_that("read", {
    for (type in c("nifti", "txt")) {
        
        vcat(T, "For %s", type)
        
        # File Setup
        nifti_files <- get_filelist(type)
        mask_files <- get_masklist(type)
        grp_mask_file <- "data/test_mask_grp.nii.gz"
        
        # Data Setup
        inlist <- load_funcs.prepare(nifti_files)
        inlist <- load_funcs.mask(inlist, automask=TRUE, subject.masks=mask_files, 
                                  group.mask=grp_mask_file)
        
        # Reference
        ref.inlist <- inlist
        inds <- which(inlist$mask)
        ref.inlist$funcs <- lapply(inlist$files, function(f) {
            if (type == "nifti") {
                func <- inlist$reader(f, voxs=inds)
            } else {
                func <- inlist$reader(f)
                func <- deepcopy(func, cols=inds)
            }
        })
        ref <- lapply(ref.inlist$funcs, as.matrix)
        
        # Comparison
        comp.inlist <- load_funcs.read(inlist)
        comp <- lapply(comp.inlist$funcs, as.matrix)
        
        expect_that(ref, equals(comp))
    }
})






context("load_funcs.scale")

test_that("scale", {
    for (type in c("nifti", "txt")) {
        
        vcat(T, "For %s", type)
        
        # File Setup
        nifti_files <- get_filelist(type)
        mask_files <- get_masklist(type)
        grp_mask_file <- "data/test_mask_grp.nii.gz"
        
        # Data Setup
        inlist <- load_funcs.prepare(nifti_files)
        inlist <- load_funcs.mask(inlist, automask=TRUE, subject.masks=mask_files, 
                                  group.mask=grp_mask_file)
        inlist <- load_funcs.read(inlist)
        
        # Reference
        ref.inlist <- inlist
        ref.inlist$funcs <- lapply(ref.inlist$funcs, scale)
        ref <- lapply(ref.inlist$funcs, as.matrix)
        
        # Comparison
        comp.inlist <- load_funcs.scale(inlist)
        comp <- lapply(comp.inlist$funcs, as.matrix)
        
        expect_that(ref, equals(comp))
    }
})

test_that("scale works when copying func", {
    for (type in c("nifti", "txt")) {
        
        vcat(T, "For %s", type)
        
        # File Setup
        nifti_files <- get_filelist(type)
        mask_files <- get_masklist(type)
        grp_mask_file <- "data/test_mask_grp.nii.gz"
        
        # Data Setup
        inlist <- load_funcs.prepare(nifti_files)
        inlist <- load_funcs.mask(inlist, automask=TRUE, subject.masks=mask_files, 
                                  group.mask=grp_mask_file)
        inlist <- load_funcs.read(inlist)
        
        # Reference
        ref.inlist <- inlist
        ref.inlist$funcs_nonscaled <- ref.inlist$funcs
        ref.inlist$funcs <- lapply(ref.inlist$funcs, scale, to.copy=TRUE)
        ref1 <- lapply(ref.inlist$funcs, as.matrix)
        ref2 <- lapply(ref.inlist$funcs_nonscaled, as.matrix)
        
        # Comparison
        comp.inlist <- load_funcs.scale(inlist, to.copy=TRUE)
        comp1 <- lapply(comp.inlist$funcs, as.matrix)
        comp2 <- lapply(comp.inlist$funcs_nonscaled, as.matrix)
        
        expect_that(ref1, equals(comp1))
        expect_that(ref2, equals(comp2))
    }
})




context("load_funcs")

test_that("the big guy/gal", {
    for (type in c("nifti", "txt")) {
        
        vcat(T, "For %s", type)
        
        # File Setup
        nifti_files <- get_filelist(type)
        mask_files <- get_masklist(type)
        grp_mask_file <- "data/test_mask_grp.nii.gz"
        
        # Reference
        ref.inlist <- load_funcs.prepare(nifti_files)
        ref.inlist <- load_funcs.mask(ref.inlist, automask=TRUE, 
                                      subject.masks=mask_files, 
                                      group.mask=grp_mask_file)
        ref.inlist <- load_funcs.read(ref.inlist)
        ref.inlist <- load_funcs.scale(ref.inlist, TRUE, FALSE, FALSE)
        ref <- lapply(ref.inlist$funcs, as.matrix)
        
        # Comparison
        comp.inlist <- load_funcs(nifti_files, automask=TRUE, 
                                  subject.masks=mask_files, 
                                  group.mask=grp_mask_file)
        comp <- lapply(comp.inlist$funcs, as.matrix)
        
        expect_that(ref.inlist$files, equals(comp.inlist$files))
        expect_that(ref.inlist$ftype, equals(comp.inlist$ftype))
        expect_that(ref.inlist$mask, equals(comp.inlist$mask))
        expect_that(ref, equals(comp))
        expect_that(ref.inlist$reader, equals(comp.inlist$reader))
    }
})

