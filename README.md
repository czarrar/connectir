Connectome-Wide Association 
=========

Advanced functional connectivity functions and scripts including ones for connectome-wide association studies. This package offers a fast and memory efficient approach to analyzing big brain connectivity data. It also provides a framework for analyzing how variations in an n x n connectivity matrix can explain behavior or psychiatric diagnosis (where n = a node or brain region). The approach first examines the most important nodes in the brain associated with behavior. It then examines what connections in those relevant nodes are driving the association with behavior. This hierarchical approach makes the connectome-wide association results more interpretable.

> If you are interested in a permutation-based cluster correction of data, please see this gist: https://gist.github.com/czarrar/0e620cf64a975380cdc2.

> Please note that this package is no longer actively maintained. You can get it to run using an older version of R, Rcpp, and other packages. However, I have ported the R code into Python. Please use the Python version of the package that is a part of [CPAC](https://github.com/FCP-INDI/C-PAC).

# Overview

Connectome-wide association studies (CWAS) remain limited by statistical approaches that are computationally intensive, depend on a priori hypotheses, or require stringent correction for multiple comparisons. To address these issues, we proposed a data-driven framework (Shehzad et al., 2014) that provides a comprehensive, voxel-wise survey of brain-behavior relationships using multivariate distance matrix regression (MDMR; Anderson 2001). The approach identifies voxels whose whole-brain connectivity patterns vary significantly with a phenotypic variable. Here, we introduce the R software package `connectir` for applying CWAS on resting-state fMRI data.

# Installation

## Quick Approach

1. Install R and optionally Rstudio.
2. Install the relevant packages within R including connectir using my script [connectir_install.R](https://github.com/czarrar/Rinstall/blob/master/connectir_install.R).

## Details and Troubleshooting

### Parallel Matrix Algebra Operations

There are two ways to parallelize the analyses. One approach is to divide your workflow into smaller chunks and run those separately (like separate processes). This comes with the R packages installed with connectir_install.R. Another approach is to run each matrix algebra operation (e.g., dot product) in parallel, which we go into detail in this section. Below I detail different linear algebra libraries and linking them to R. Note this section is still under development.

#### Intel MKL

If you have Windows, Ubuntu, or RedHat/Centos, you can install Revolution R. This is a version of R compiled with Intel MKL by the company Revolution Analytics available free for academic use. You can get it from here.

Another option is to compile and install R linked with Intel MKL on your own. Here is a good and quick tutorial.

#### OpenBlas

You can also install R via my own script that links R with a parallel matrix algebra library called openblas. This script is in the Rinstall repo and is called install.py.

Another option for linux is to download repositories. A good/quick tutorial can be found here.

### Installing Connectir and Other R Packages

After R is setup, there are several packages within R that need to be installed. To do this, please run connectir_install.R. After downloading (or copying and pasting) this script to your machine, you can run it with Rscript [connectir_install.R](https://github.com/czarrar/Rinstall/blob/master/connectir_install.R). On certain linux systems, you need to ensure you have libcurl and libxml installed.


# Tutorial

Here we give a vanilla run of CWAS-MDMR and further details can be found on the [wiki](https://github.com/czarrar/connectir/wiki). I also go through these steps in our [resting-state conference poster](https://www.dropbox.com/s/5thqpxk7a9lueis/zarrar_rsn_poster_2014_v1.pdf?dl=0).

## Subject Distances

```bash
connectir_subdist.R \ 
-i functional_path_list.txt \ 
--automask1 \
 --brainmask1 standard_grey_matter.nii.gz \ 
--bg standard_brain_4mm.nii.gz \ 
--memlimit 20 -c 3 -t 4 \
 subject_distances_outdir
 ```
 
 * -i: List of your input functional images. Can be nifti (nii or nii.gz) containing voxelwise time-series or text files containing region/parcellation time-series (columns=regions and rows=time-points).
* --automask1: Will generate the group mask containing only voxels that have non-zero values (i.e., variance) across all participants.
* --brainmask1: An additional prior mask. We tend to use a 25% probability grey matter mask in MNI152 standard space. You can find these on the CPAC website.
* --bg: This is used to determine writing of output voxelwise files and also in the future will be used to generate images of the results. Since my data here is assumed to be voxelwise in 4mm space, the standard reference image is also in 4mm space.
* --memlimit: The memory (RAM) limit of the processes in GB. Here it is set to 20GB.
* -c: Number of parallel jobs/forks to run in parallel.
* -t: Number of parallel linear algebra operations.
* Finally the last argument gives the full path to the output directory.

## Multivariate Distance Matrix Regression (MDMR)

For the options that are the same as before (--memlimit, -c, -t), I will not repeat the description here.

```bash
connectir_mdmr.R \
 -i subject_distances_outdir \
 --formula FSIQ + Age + Sex + meanFD \
 --model model_evs.csv \ 
--factors2perm FSIQ \ 
--memlimit 8 -c 3 -t 4 \
 --save-perms --ignoreprocerror \
 iq_outdir.mdmr
```

* -i: Input path to your subject distances directory
* --formula:: We use R's formula syntax. Each variable represents a column in your model file (i.e., a factor or explanatory variable). The + indicates to combine the different variables in one model. If you want to do an interaction you can use `*`, which will generate the main effects and interaction. To only look at the interaction (for instance between Age and Sex), you would use : as in Age:Sex.
* --model: This is your model file containing all your explanatory variables and covariates in each column. Each row would correspond to each subject/scan in the same order as the input functional file list in connectir_subdist.R. So ensure that this model file as the rows (scans/subjects) correspond to the rows in -i file list of connectir_subdist.R.
* --factors2perm: This indicates the variables in your formula to permute and calculate significance estimates. The other variables will be treated as covariates. In this case, we were interested in full-scale IQ (FSIQ).
* --saveperms: This outputs all the permuted F-statistics into a file. This is needed for now to later calculate (in a separate script), permuted cluster correction.
* --ignoreprocerror: Sometimes the code tries to estimate the maximum number of cores on your machine and does this wrong. This option ignores the error that's thrown because of this potential issue.
* The final argument is the output directory for the MDMR script. This should just be the directory name (not the path). The path is set to the subject distances directory given with the -i option.

# Publications

Please refer to the [paper in NeuroImage](http://www.sciencedirect.com/science/article/pii/S1053811914001232)
