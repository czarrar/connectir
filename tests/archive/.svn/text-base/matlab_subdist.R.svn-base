library(connectir)
setwd("/mnt/nfs/class/exampledata/fmri/adhd200/command")

# settings
opts <- list(
    forks = 1,
    threads = 12, 
    verbose = TRUE, 
    memlimit = 40, 
    blocksize = 0, 
    brainmask = "/mnt/nfs/class/exampledata/fmri/adhd200/data/niak/rois_3000_mask.txt", 
    infuncs = as.character(read.table("tseries_fnames.txt")[,1]), 
    ztransform = TRUE,
    design_file = "xxx"
)

# set parallel processing
set_parallel_procs(opts$forks, opts$threads, opts$verbose)
parallel <- ifelse(opts$forks>1, TRUE, FALSE)

# read design matrix
X <- as.matrix(read.table("design.txt", header=T))
if (length(opts$infuncs) != nrow(X))
    stop("size mismatch")
design <- as.big.matrix(X, shared=parallel, type="double")
r <- qlm_rank(design)
if (r < ncol(design))
    stop("design matrix is rank deficient")

# read mask
mask <- as.logical(read.table(opts$brainmask)[,1])
nvoxs <- sum(mask)

# load the data and mask
library(R.matlab)
reader <- gen_big_reader("matlab")
fnames <- opts$infuncs

#fnames <- opts$infuncs[1:10]
#design <- as.big.matrix(design[1:10,])
#design <- as.big.matrix(design[,c(1,10:12)])

funcs <- load_and_mask_func_data2(fnames, reader, mask, verbose=opts$verbose, parallel=parallel)
tpts <- sapply(funcs, nrow)

# calculate memlimit
nsubs <- length(funcs)
opts <- get_subdist_memlimit(opts, nsubs, nvoxs, tpts)

#opts$blocksize <- 1000

# create subject distance matrix
inds <- 1:nvoxs
outdir <- "/mnt/nfs/class/exampledata/fmri/adhd200/analysis/niak_all_all_run1_rois3000_895subs.subdist"

# cor maps => regress out some covariates => subject distances
subdist <- big.matrix(nsubs^2, nvoxs, type="double", shared=parallel)
compute_subdist2(funcs, subdist, inds, opts$blocksize, opts$ztransform, 
                  verbose=opts$verbose, design_matrix=design)
tmp <- deepcopy(subdist, backingpath=outdir, backingfile="subdist_residuals.bin", descriptorfile="subdist_residuals.desc")
rm(subdist); rm(tmp); gc(FALSE, TRUE)

# new 
subdist <- big.matrix(nsubs^2, nvoxs, type="double", shared=parallel)
compute_subdist2(funcs, subdist, inds, opts$blocksize, opts$ztransform, 
                  verbose=opts$verbose)
tmp <- deepcopy(subdist, backingpath=outdir, backingfile="subdists.bin", descriptorfile="subdists.desc")
rm(subdist); rm(tmp); gc(FALSE, TRUE)


compute_subdist2(funcs, subdist, inds, opts$blocksize, opts$ztransform, 
                 verbose=opts$verbose)



# save (?gower?)


n <- 1000
x <- matrix(rnorm(n*n), n, n)
y <- matrix(rnorm(n*n), n, n)
system.time(x %*% y)

bx <- as.big.matrix(x)
by <- as.big.matrix(y)
bz <- big.matrix(n, n)
system.time(dgemm(C=bz, A=bx, B=by))

require(bigmemory)
A = matrix(rnorm(100), 10)
B = matrix(rnorm(100), 10)
C = matrix(0, 10, 10)
print(A %*% B)
dgemm(C=C, A=A, B=B)
print(C)
A = as.big.matrix(A, type='double')
C = as.big.matrix(C, type='double')
dgemm(C=C, A=A, B=B)
print(head(C, 10))
