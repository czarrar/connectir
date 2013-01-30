library(connectir)
bpath <- "/mnt/nfs/class/exampledata/fmri/adhd200/analysis/niak_all_all_run1_rois3000_895subs.subdist"
setwd(bpath)

# General Function(s)
printf <- function(msg, ..., newline=TRUE) {
    if (opts$verbose) {
        cat(sprintf(msg, ...))
        if (newline) cat("\n")
    }
}

Sys.setenv(MKL_NUM_THREADS=10)



# 1. Load stuff

save_pvals <- c()
x1 <- attach.big.matrix("subdist_residuals.desc")
x2 <- attach.big.matrix("subdists.desc")

demog <- read.csv("demog_run1_all_included.txt", sep="\t")
demog.train <- subset(demog, type=="training")
which.train <- which(demog$type=="training")
dtf <- data.frame(demog.train[,c('DX', 'Full4.IQ', 'Age', 'Gender')])
dtf$DX <- factor(as.character(dtf$DX))

dtf_adhdc_tdc <- dtf[dtf$DX=='TD'|dtf$DX=='ADHD-C',]
dtf_adhdc_tdc$DX <- factor(as.character(dtf_adhdc_tdc$DX))
w_adhdc_tdc <- which.train[dtf$DX=='TD'|dtf$DX=='ADHD-C']

dtf_adhdi_tdc <- dtf[dtf$DX=='TD'|dtf$DX=='ADHD-I',]
dtf_adhdi_tdc$DX <- factor(as.character(dtf_adhdi_tdc$DX))
w_adhdi_tdc <- which.train[dtf$DX=='TD'|dtf$DX=='ADHD-I']

dtf_adhdc_adhdi <- dtf[dtf$DX=='ADHD-C'|dtf$DX=='ADHD-I',]
dtf_adhdc_adhdi$DX <- factor(as.character(dtf_adhdc_adhdi$DX))
w_adhdc_adhdi <- which.train[dtf$DX=='ADHD-C'|dtf$DX=='ADHD-I']

f1 <- . ~ DX
f2 <- . ~ DX + Full4.IQ + Age + Gender

which_rois <- which(as.logical(read.table("rois_3000_mask.txt")[,1]))
mask <- read.mask("rois_3000_brainmask.nii.gz")
hdr <- read.nifti.header("rois_3000_brainmask.nii.gz")
rois <- read.mask("rois_3000.nii.gz", NULL)
rois <- rois[mask]


# 2. ADHD-C vs. TDC

model <- dtf_adhdc_tdc
xdist <- slice.subdist(x1, subs=w_adhdc_tdc, shared=FALSE)  # will create a copy
x1 <- free.memory(x1, bpath)
invisible(gc(FALSE, TRUE))

gdist <- gower.subdist(xdist, shared=FALSE)
rm(xdist); invisible(gc(FALSE, TRUE))


## BLOCK SIZE

opts <- list()
opts$verbose <- TRUE
opts$permutations <- 5000
opts$blocksize <- 0
opts$memlimit <- 10
opts$factors2perm <- "DX"

# functions
n2gb <- function(x) x*8/1024^3
gb2n <- function(x) x/8*1024^3

# general info
nsubs <- sqrt(nrow(gdist))
nvoxs <- ncol(gdist)
if (nrow(model) != nsubs)
    stop("Number of subjects in distance matrix does not match number of subjects in model")

# RAM for distance matrix
mem_used4dmat <- n2gb(nsubs^2 * nvoxs)

# RAM for Fperms matrix
mem_used4fperms <- n2gb(opts$permutations * nvoxs)

if (opts$blocksize==0) {
    printf("...autosetting blocksize to -> ", newline=F)
  
    # minimum amount of RAM needed (assume blocksize of 2)
    min_mem_needed <- n2gb(
        opts$permutations * 2 * 2   # for temporary matrices
        + nsubs^2 * 2
    ) * getDoParWorkers() + mem_used4dmat + mem_used4fperms
  
    # limit in RAM use
    mem_limit <- as.numeric(opts$memlimit)
    if (mem_limit < min_mem_needed)
        stop(sprintf("You require at least %.2f GB of memory but are limited to %i GB. Please set the --memlimit option to a higher number to continue.", min_mem_needed, mem_limit))
  
    # amount of RAM for mdmr
    mem_used4mdmr <- mem_limit - mem_used4dmat - mem_used4fperms
  
    # block size
    opts$blocksize <- floor(gb2n(
        mem_used4mdmr/((4*opts$permutations+2*nsubs^2)*getDoParWorkers())
    ))
    printf("%i", opts$blocksize)
  
    # clear variables
    rm(min_mem_needed, mem_limit, mem_used4mdmr)
  
}

# temporary RAM needed while computing MDMR
## error variance and explained variance
mem_used4temp <- n2gb(opts$permutations * opts$blocksize * 2 * getDoParWorkers())

# copy of subdist
mem_used4chunk <- n2gb(nsubs^2 * opts$blocksize * getDoParWorkers())

# total ram
mem_used <- mem_used4dmat + mem_used4fperms + mem_used4temp + mem_used4chunk

printf("...will be using %.2f GB of RAM", mem_used)
rm(nsubs, nvoxs, mem_used4dmat, mem_used4fperms, mem_used4temp, mem_used4chunk, mem_used)
invisible(gc(FALSE))


# Compute MDMR
printf("05. Computing MDMR")
start.time <- Sys.time()

# NOTE: need to paste new version!

res.mdmr2 <- mdmr(gdist, f1, model, nperms=opts$permutations, block.size=opts$blocksize, verbose=opts$verbose, factors2perm=opts$factors2perm, strata=NULL, shared=TRUE)
rm(gdist)
invisible(gc(FALSE, TRUE))

end.time <- Sys.time()
printf("MDMR is done! It took: %.2f minutes\n", as.numeric(end.time-start.time, units="mins"))

pvals <- res.mdmr2$pvals[,1]
nii_pvals <- rois*0
for (i in 1:length(which_rois))
    nii_pvals[rois==which_rois[i]] = 1-pvals[i]
write.nifti(nii_pvals, hdr, mask, outfile="mdmr_residuals_adhdc_vs_tdc_1p.nii.gz")

save_pvals <- cbind(save_pvals, pvals)

rm(res.mdmr2)


# 3. ADHD-I vs. TDC

model <- dtf_adhdi_tdc
xdist <- slice.subdist(x1, subs=w_adhdi_tdc, shared=FALSE)  # will create a copy
x1 <- free.memory(x1, bpath)
invisible(gc(FALSE, TRUE))

gdist <- gower.subdist(xdist, shared=FALSE)
rm(xdist); invisible(gc(FALSE, TRUE))

# Compute MDMR
printf("05. Computing MDMR")
start.time <- Sys.time()

# NOTE: need to paste new version!

res.mdmr2 <- mdmr(gdist, f1, model, nperms=opts$permutations, block.size=opts$blocksize, verbose=opts$verbose, factors2perm=opts$factors2perm, strata=NULL, shared=TRUE)
rm(gdist)
invisible(gc(FALSE, TRUE))

end.time <- Sys.time()
printf("MDMR is done! It took: %.2f minutes\n", as.numeric(end.time-start.time, units="mins"))

pvals <- res.mdmr2$pvals[,1]
nii_pvals <- rois*0
for (i in 1:length(which_rois))
    nii_pvals[rois==which_rois[i]] = 1-pvals[i]
write.nifti(nii_pvals, hdr, mask, outfile="mdmr_residuals_adhdi_vs_tdc_1p.nii.gz")

save_pvals <- cbind(save_pvals, pvals)

rm(res.mdmr2)


# 3. ADHD-C vs. ADHD-I

model <- dtf_adhdc_adhdi
xdist <- slice.subdist(x1, subs=w_adhdc_adhdi, shared=FALSE)  # will create a copy
x1 <- free.memory(x1, bpath)
invisible(gc(FALSE, TRUE))

gdist <- gower.subdist(xdist, shared=FALSE)
rm(xdist); invisible(gc(FALSE, TRUE))

# Compute MDMR
printf("05. Computing MDMR")
start.time <- Sys.time()

# NOTE: need to paste new version!

res.mdmr2 <- mdmr(gdist, f1, model, nperms=opts$permutations, block.size=opts$blocksize, verbose=opts$verbose, factors2perm=opts$factors2perm, strata=NULL, shared=TRUE)
rm(gdist)
invisible(gc(FALSE, TRUE))

end.time <- Sys.time()
printf("MDMR is done! It took: %.2f minutes\n", as.numeric(end.time-start.time, units="mins"))

pvals <- res.mdmr2$pvals[,1]
nii_pvals <- rois*0
for (i in 1:length(which_rois))
    nii_pvals[rois==which_rois[i]] = 1-pvals[i]
write.nifti(nii_pvals, hdr, mask, outfile="mdmr_residuals_adhdc_vs_adhdi_1p.nii.gz")

save_pvals <- cbind(save_pvals, pvals)

rm(res.mdmr2)

