#' Correct MDMR results for multiple comparisons
#' based on maximum cluster size at each permutation.
#' Assumes that you have proper mdmr output directory.
#' 
#' Will save the significant clusters, p-values masked by the significant 
#' clusters, and Z-scores masked by the significant clusters.
#'
#' @author Zarrar Shehzad
#' @param mdmr.dir MDMR output directory
#' @param factor name of factor or predictor variable
#' @param ... see \code{\link{clust_mdmr.correct}}
#' @return vector of non-zero voxels indicating significant clusters
#' @seealso \code{\link{clust_mdmr.correct}}
clust_mdmr.correct_wrapper <- function(mdmr.dir, factor, ...)
{
    # Set paths
    sdist.dir <- dirname(mdmr.dir)
    maskfile <- file.path(sdist.dir, "mask.nii.gz")
    fstatsfile <- file.path(mdmr.dir, sprintf("fperms_%s.desc", factor))
    pvalsfile <- file.path(mdmr.dir, sprintf("pvals_%s.nii.gz", factor))
    
    # Check paths
    if (!file.exists(mdmr.dir))
        vstop("MDMR directory '%s' cannot be found", mdmr.dir)
    if (!file.exists(maskfile))
        vstop("Mask '%s' cannot be found", maskfile)
    if (!file.exists(fstatsfile))
        vstop("F-statistics '%s' cannot be found", fstatsfile)
    if (!file.exists(pvalsfile))
        vstop("P-values '%s' cannot be found", pvalsfile)
    
    # Read in mask, header, fstats, and p-values
    mask <- read.mask(maskfile)
    hdr <- read.nifti.header(maskfile)
    fstats <- attach.big.matrix(fstatsfile)
    pvals <- read.nifti.image(pvalsfile)[mask]
        
    # Cluster correct
    clust <- clust_mdmr.correct(mask, hdr, fstats, mdmr.dir, ...)
    
    # Save significant clusters
    ofile <- file.path(mdmr.dir, sprintf("clust_%s.nii.gz", factor))
    write.nifti(clust, hdr, mask, odt="int", outfile=ofile)
    
    # Save p-values with only significant clusters
    ofile <- file.path(mdmr.dir, sprintf("clust_pvals_%s.nii.gz", factor))
    clust_pvals <- pvals * (clust>0)
    write.nifti(clust_pvals, hdr, mask, odt="float", outfile=ofile)
    
    # Save Z-scores for significant clusters
    ofile <- file.path(mdmr.dir, sprintf("clust_zstats_%s.nii.gz", factor))
    clust_zstats <- qt(pvals, Inf) * (clust>0)
    write.nifti(clust_zstats, hdr, mask, odt="float", outfile=ofile)
    
    return(clust)
}

#' Correct MDMR results for multiple comparisons
#' based on maximum cluster size at each permutation
#'
#' @author Zarrar Shehzad
#' @param mask boolean values indicating non-zero voxels
#' @param hdr list header usually from brain mask image
#' @param fstats big.matrix of nperms x nvoxs
#' @param fpath file path (if fstats is file.backed)
#' @param vox.thresh voxel-level threshold (p-value)
#' @param clust.thresh cluster-level threshold (p-value)
#' @param clust.type type of cluster threshold (mass or size)
#' @param verbose boolean
#' @param parallel boolean
#' @return vector of non-zero voxels indicating significant clusters
clust_mdmr.correct <- function(mask, hdr, fstats, fpath=NULL, 
                            vox.thresh=0.05, clust.thresh=0.05, 
                            clust.type=c("mass", "size"),  
                            verbose=TRUE, parallel=FALSE)
{
    vcat(verbose, "Multiple comparisons correction with cluster %s", 
            clust.type)
    
    # Get p-values for each permutation
    vcat(verbose, "\t...getting p-values")
    pvals <- clust_mdmr.get_pvals_for_each_perm(fstats, fpath, verbose, 
                                                parallel)
    
    # Get clust size/masses
    vcat(verbose, "\t...getting cluster sizes")
    clust_vals <- clust_mdmr.values(pvals, mask, hdr$dim, vox.thresh, 
                                    verbose, parallel)
    
    # Get cluster p-values
    vcat(verbose, "\t...getting cluster p-values")
    ref_details <- cluster.table(1-pvals[1,], 1-vox.thresh, hdr$dim, mask)
    ref_details <- clust_mdmr.pvalues(clust_vals, ref_details, clust.type)
    
    # Get clusters on brain image
    vcat(verbose, "\t...getting clusters")
    clust <- clust_mdmr.clusterize(ref_details, mask, clust.thresh)
    
    return(clust)
}

#' Extract a column from a big matrix
#'
#' @author Zarrar Shehzad
#' @param x big.matrix
#' @param i column index
#' @param bp backing-path if x is file.backed
#' @return desired column as a vector
bm.get_col <- function(x, i, bp=NULL) {
    col <- matrix(0, nrow(x), 1)
    subx <- sub.big.matrix(x, firstCol=i, lastCol=i, backingpath=bp)
    dcopy(X=subx, Y=col)
    as.vector(col)
}

#' Save a vector as a column in a big matrix
#'
#' @author Zarrar Shehzad
#' @param x column vector
#' @param y big.matrix
#' @param i column index in big.matrix
#' @param bp backing-path if y is file.backed
#' @return big matrix y
bm.save_col <- function(x, y, i, bp=NULL) {
    x <- matrix(x, length(x), 1)
    suby <- sub.big.matrix(y, firstCol=i, lastCol=i, backingpath=bp)
    deepcopy(x=x, y=suby)
}

#' Get p-values for each permutation
#'
#' @author Zarrar Shehzad
#' @param fstats big.matrix of F-statistics, nperms x nvoxs
#' @param fpath if fstats is file.backed, path to fstats
#' @param verbose TRUE or FALSE
#' @param parallel TRUE or FALSE (must also set nforks)
#' @return big.matrix of nperms x nvoxs containing p-values
clust_mdmr.get_pvals_for_each_perm <- function(fstats, fpath=NULL,
                                            verbose=T, parallel=F)
{
    if (is.filebacked(fstats) && is.null(fpath))
        vstop("fstats is file-backed, must supply fpath")
    if (!is.null(fpath) && !file.exists(fpath))
        vstop("could not find fpath: %s", fpath)
    
    progress <- ifelse(verbose, "text", "none")
    nperms <- nrow(fstats)
    nvoxs <- ncol(fstats)
    
    ps.mat <- big.matrix(nperms, nvoxs, type="double", shared=TRUE)
    tmp <- llply(1:nvoxs, function(i) {
        x <- bm.get_col(fstats, i, fpath)
        rx <- nperms - rank(x, ties.method="max")
        ps <- rx/nperms
        bm.save_col(ps, ps.mat, i)
        return(NULL)
    }, .progress=progress, .parallel=parallel)
    rm(tmp)
    
    return(ps.mat)
}

#' Get maximum cluster mass, size, and mass/size
#' for each permutation's voxelwise p-value map
#'
#' @author Zarrar Shehzad
#' @param ps.mat big.matrix of p-values, nperms x nvoxs
#' @param mask boolean vector of non-zero voxels
#' @param voxdim vector with the size of each dimension for mask 
#'        and each ps.mat column
#' @param vox.thresh voxel-level threshold (p-value)
#' @param verbose TRUE or FALSE
#' @param parallel TRUE or FALSE (must also set nforks)
#' @return array of nperms x 3 containing max cluster values
#'         for each permutation
clust_mdmr.values <- function(ps.mat, mask, voxdim, 
                                vox.thresh=0.05, 
                                verbose=T, parallel=F)
{
    progress <- ifelse(verbose, "text", "none")
    nperms <- nrow(ps.mat)
    nvoxs <- ncol(ps.mat)
    
    clust_values <- laply(1:nperms, function(i) {
        ct <- cluster.table(1-ps.mat[i,], 1-vox.thresh, voxdim, mask)
        c(size=ct$max.size, mass=ct$max.mass, rel=ct$max.rel)
    }, .progress=progress, .parallel=parallel)
    
    return(clust_values)
}

#' Gets the p-values for each cluster
#' based on the cluster's mass or size
#'
#' @author Zarrar Shehzad
#' @param clust_vals is a matrix with size, mass, & mass/size values 
#'        and is nperms x 3
#' @param ref_details is a list resulting from cluster.table
#'        and it has various cluster details
#' @param clust.type mass or size
#' @return ref_details but with cluster p-values added to the list
clust_mdmr.pvalues <- function(clust_vals, ref_details, 
                                clust.type=c("mass", "size"))
{
    nperms <- nrow(clust_vals)
    
    if (clust.type == "size") {
        clust.pvals <- sapply(ref_details$size, function(x) 
                                sum(clust_vals[,1]>=x))/nperms
    } else if (clust.type == "mass") {
        clust.pvals <- sapply(ref_details$mass, function(x) 
                                sum(clust_vals[,2]>=x))/nperms
    } else {
        vstop("unrecognized clust.type: %s", clust.type)
    }
    
    ref_details$clust.pvals <- clust.pvals
    
    return(ref_details)
}

#' Generate's brain image with significant clusters only
#' 
#' @author Zarrar Shehzad
#' @param ref_details is a list resulting from cluster.table
#'        and clust_mdmr.pvalues
#' @param mask brain mask
#' @param clust.thresh p-value for cluster threshold
#' @return brain image with significant clusters (only voxels in mask)
clust_mdmr.clusterize <- function(ref_details, mask, clust.thresh=0.05)
{
    # Get clusters
    clust <- ref_details$clust[mask]
    
    # Mark significant clusters
    w.clusts <- which(rev(ref_details$clust.pvals<clust.thresh))
    
    # Create new image with significant clusters as numbers
    clust.numbered <- clust*0    # empty vector
    for (i in 1:length(w.clusts)) clust.numbered[clust==w.clusts[i]] <- i
    
    return(clust.numbered)
}

# obj list: nfactors, fpath, fstats, Pmat
clust_mdmr <- function(obj, maskfile, vox.thresh=0.05, clust.thresh=0.05, 
                        clust.type=c("mass", "size"), verbose=TRUE, parallel=FALSE)
{
    vcat(verbose, "Multiple comparisons correction with cluster %s", clust.type)
    if (!(clust.type %in% c("mass", "size")))
        vstop("unrecognized clust.type %s", clust.type)
    progress <- ifelse(verbose, "text", "none")
    nvoxs <- nrow(obj$Pmat)
    nfactors <- obj$nfactors
    fpath <- obj$fpath
    mask <- read.mask(maskfile)
    hdr <- read.nifti.header(maskfile)
    
    get_col <- function(x, i, bp=NULL) {
        col <- matrix(0, nrow(x), 1)
        subx <- sub.big.matrix(x, firstCol=i, lastCol=i, backingpath=bp)
        dcopy(X=subx, Y=col)
        as.vector(col)
    }
    
    save_col <- function(x, y, i, bp=NULL) {
        x <- matrix(x, length(x), 1)
        suby <- sub.big.matrix(y, firstCol=i, lastCol=i, backingpath=bp)
        deepcopy(x=x, y=suby)
    }
    
    Cmat <- big.matrix(nvoxs, nfactors, type="double", shared=TRUE)
    
    refs <- c()
    for (fi in 1:nfactors) {
        vcat(verbose, "\tFactor #%i",fi )
        # Get inputs
        pvals <- obj$Pmat[,fi]
        fs.mat <- obj$fstats[[fi]]
        nperms <- nrow(fs.mat)
        nvoxs <- ncol(fs.mat)
        
        # Get p-values for each permutation
        vcat(verbose, "\t...getting p-values")
        ps.mat <- big.matrix(nperms, nvoxs, type="double", shared=parallel)
        tmp <- llply(1:nvoxs, function(i) {
            x <- get_col(fs.mat, i, fpath)
            rx <- nperms - rank(x, ties.method="max")
            ps <- rx/nperms
            save_col(ps, ps.mat, i)
            return(NULL)
        }, .progress=progress, .parallel=parallel)
        rm(tmp)
        
        # Get cluster sizes/masses
        vcat(verbose, "\t...getting cluster sizes")
        ref <- cluster.table(1-pvals, 1-vox.thresh, hdr$dim, mask)
        comps <- laply(1:nperms, function(i) {
            ct <- cluster.table(1-ps.mat[i,], 1-vox.thresh, hdr$dim, mask)
            c(size=ct$max.size, mass=ct$max.mass, rel=ct$max.rel)
        }, .progress="text", .parallel=parallel)
        
        # P-values for cluster size/mass
        vcat(verbose, "\t...getting cluster p-values and clusters")
        if (clust.type == "size") {
            clust.pvals <- sapply(ref$size, function(x) sum(comps[,1]>=x))/nperms
        } else {
            clust.pvals <- sapply(ref$mass, function(x) sum(comps[,2]>=x))/nperms
        }
        ref$clust.pvals <- clust.pvals
        refs[[fi]] <- ref
        
        # Clusters
        clust <- ref$clust[mask]
        w.clusts <- which(rev(clust.pvals<clust.thresh))
        clust.new <- clust*0
        for (i in 1:length(w.clusts)) clust.new[clust==w.clusts[i]] <- i
        Cmat[,fi] <- clust.new
        
        # Clear f permutations from memory
        if (is.filebacked(fs.mat)) {
            obj$fstats[[fi]] <- free.memory(fs.mat, fpath)
        }
    }
    
    return(list(mat=Cmat, ref=refs))
}
