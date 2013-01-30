# This assumes that connectir_subdist has been run once (annoying I know)

load("options.rda")

printf <- function(msg, ..., newline=TRUE) {
    if (opts$verbose) {
        cat(sprintf(msg, ...))
        if (newline) cat("\n")
    }
}

#####
# MASKS
#####

brainmask <- read.mask("input_masks/brainmask.nii.gz")
seedmask <- read.mask("input_masks/seedmask.nii.gz")


#####
# FUNCTIONAL DATA
#####

infiles <- as.character(opts$infiles)
funclist <- load_and_mask_func_data(infiles, brainmask, shared=FALSE, type="double")


#####
# Distance Matrix
#####

nsubs <- length(funclist)
nvoxs <- sum(brainmask)
subdist <- big.matrix(nsubs^2, nvoxs, type="double")


#####
# Test
#####

inds <- 1:10
cormaps_list <- vbca_batch(funclist, 1:10, ztransform=ztransform, shared=FALSE)
subdist_CHUNK <- sub.big.matrix(subdist, firstCol=1, lastCol=10)
tmp <- compute_subdist_worker(cormaps_list, inds, subdist_CHUNK)
