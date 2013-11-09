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

load_inlist <- function(inlist, scale=FALSE, parallel=TRUE) {
    funclist <- load_and_mask_func_data2(inlist$files, inlist$reader, mask=inlist$mask, 
                                         verbose=TRUE, scale=scale, 
                                         type="double", shared=parallel)
    return(funclist)
}
