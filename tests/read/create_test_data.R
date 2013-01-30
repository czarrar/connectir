# This isn't really for testing but code I used to create some of the test functionals

library(niftir)

funcfile <- system.file("data/test_func_inds.nii.gz", package="niftir")
hdr <- read.nifti.header(funcfile)
hdr$dim <- hdr$dim[1:3]; hdr$pixdim <- hdr$pixdim[1:3]
mask <- rep(1, prod(hdr$dim))
inds1 <- sample(1:length(mask))[1:2]
inds2 <- sample(1:length(mask))[1:2]
for (i in 1:9) {
    if (i %% 2) {
        inds <- inds1
    } else {
        inds <- inds2
    }
    
    subj.mask <- mask
    subj.mask[inds] <- 0
    
    ofile <- sprintf("data/test_mask%02i.nii.gz", i)
    write.nifti(subj.mask, hdr, outfile=ofile, overwrite=T)
    
    ofile <- sprintf("data/test_mask%02i.txt", i)
    write.table(subj.mask, file=ofile, row.names=F, col.names=F, quote=F)
    
    ifile <- system.file("data/test_func_inds.nii.gz", package="niftir")
    ofile <- sprintf("data/test_func%02i.nii.gz", i)
    bm <- read.big.nifti(ifile, nifti4d=TRUE, type="double")
    bm[,] <- rnorm(sum(bm@mask))
    bm[,inds] <- 0
    
    write.nifti(bm, outfile=ofile, odt="float", overwrite=TRUE)
    
    outfile <- sprintf("data/test_func%02i.1D", i)
    write.table(bm[,], file=outfile, row.names=F, col.names=F)
    
    infile <- sprintf("data/test_func%02i.1D", i)
    outfile <- sprintf("data/test_func%02i.txt", i)
    file.remove(outfile)
    file.copy(infile, outfile)
}

grp.mask <- mask
grp.mask[c(1,2)] <- 0
write.nifti(grp.mask, hdr, outfile="data/test_mask_grp.nii.gz", overwrite=T)
write.table(grp.mask, file="data/test_mask_grp.txt", row.names=F, col.names=F, quote=F)

