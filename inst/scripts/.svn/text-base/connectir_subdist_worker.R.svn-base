suppressPackageStartupMessages(library("optparse"))

# Make option list
option_list <- list(
    make_option(c("-i", "--infuncs1"), type="character", default=NULL, dest="infuncs1", help="File containing paths of different functional images (nifti or text files) in one column. Each node/voxel in these images will act as a 'seed'. (required)", metavar="file"),
    make_option("--infuncs2", type="character", default=NULL, dest="infuncs2", help="File containing paths of different functional images (nifti or text files) in one column. Each node/voxel in these images will act as a 'target', that is voxels from -i/--infuncs will be correlated with voxels in these images. The default is to have these images be the same as those specified in -i/--infuncs. (optional)", metavar="file"),
    make_option("--in2D1", action="store_true", default=FALSE, dest="in2d1", help="ask"), 
    make_option("--in2D2", action="store_true", default=FALSE, dest="in2d2", help="ask"), 
    make_option("--ztransform", action="store_true", default=FALSE, dest="ztransform", help="Fischer Z-Transform the correlations before calculating the distance between participants"),
    make_option("--brainmask1", type="character", default=NULL, help="When computing each whole-brain connectivity map, this mask will restrict which parts of the whole-brain are to be considered", metavar="file", dest="brainmask1"),
    make_option("--brainmask2", type="character", default=NULL, help="When computing each whole-brain connectivity map, this mask will restrict which parts of the whole-brain are to be considered only for --infuncs2", metavar="file", dest="brainmask2"), 
    make_option("--regress", type="character", default=NULL, help="A design matrix (space delimeted file where first row is a header) containing variables to regress out of each voxel's whole-brain connectivity maps before comparing distances between subjects", metavar="file"),
    make_option("--bg", type="character", default=NULL, help="Background image (e.g., MNI152 standard brain) upon which later results might be overlaid", metavar="file"), 
    make_option("--blocksize", type="integer", default=0, help="How many sets of voxels should be used in each iteration of computing the correlation (0 = auto) [default: %default]", metavar="number"),
    make_option("--memlimit", type="double", default=4, dest="memlimit", help="If blocksize is set to auto (--blocksize=0), then will set the blocksize to use a maximum of RAM specified by this option  [default: %default]", metavar="number"),
    make_option(c("-c", "--forks"), type="integer", default=1, help="Number of computer processors to use in parallel by forking the complete processing stream [default: %default]", metavar="number"),
    make_option(c("-t", "--threads"), type="integer", default=1, help="Number of computer processors to use in parallel by multi-threading matrix algebra operations [default: %default]", metavar="number"),
    make_option("--extrachecks", action="store_true", default=FALSE, help="Will do a more rigorous check of the input functionals before any calculations"),
    make_option("--sparseconn", action="store_true", default=FALSE, help="Computes  inverse covariance estimates when calculating connectivity maps"), 
    make_option("--sparsedist", action="store_true", default=FALSE, help="Computes 1 - inverse covariance estimates when calculating the distance between subject connectivity maps"), 
    make_option("--overwrite", action="store_true", default=FALSE, help="Overwrite output that already exists (default is not to overwrite already existing output)"),
    make_option("--no-link-functionals", action="store_true", default=FALSE, help="Will not create soft links to each of the functional images with the subdist directory"),
    make_option(c("-q", "--quiet"), action="store_false", dest="verbose", help="Print little output"), 
    make_option(c("-v", "--verbose"), action="store_true", default=TRUE, help="Print extra output [default]"),
    make_option(c("-d", "--debug"), action="store_true", default=FALSE, help="Like verbose but will also print more helpful error messages when --forks is >1")
)

# Make class/usage
parser <- OptionParser(usage = "%prog [options] output-directory", option_list=option_list, add_help_option=TRUE)

# Parse
parser_out <- parse_args(parser, positional_arguments = TRUE)
args <- parser_out$args
opts <- parser_out$options

if (length(args) < 1) {
    print_help(parser)
    quit(save="no", status=1)
}

saved_opts <- list(args=args, opts=opts)

tryCatch({

  # load connectir
  suppressWarnings(suppressPackageStartupMessages(library("connectir")))

  # parallel processing setup
  set_parallel_procs(opts$forks, opts$threads, opts$verbose)  
  # use foreach parallelization and shared memory?
  parallel_forks <- ifelse(opts$forks == 1, FALSE, TRUE)

  ###
  # Check/Setup Required Inputs
  ###
  
  start.time <- Sys.time()
  
  # Check required
  vcat(opts$verbose, "Checking required inputs")
  outdir <- abspath(args[1])
  if (file.exists(outdir) && !opts$overwrite)
      stop("Output directory '", outdir, "' already exists, you can use --overwrite")
  if (is.null(opts$infuncs1))
      stop("You must specify the -i/--infuncs option")
  if (!file.exists(opts$infuncs1))
      stop("The file specified by -i/--infuncs does not exist")
  if (is.null(opts$bg))
      stop("Please specify background image with --bg option")
  
  # Prepare input functional filenames
  infiles1 <- sapply(as.character(read.table(opts$infuncs1)[,1]), function(fp) {
      if (!file.exists(fp))
          stop("One of the input functionals does not exist: ", fp)
      abspath(fp)
  })
  n <- length(infiles1)
  
    
  ###
  # Check/Setup Optional Inputs
  ###

  vcat(opts$verbose, "Checking optional inputs")
  for (optname in c("brainmask1", "bg", "regress", "infuncs2", "brainmask2")) {
      arg <- opts[[optname]]
      if (!is.null(arg)) {
        if(!file.exists(arg))
            vstop("--%s file '%s' does not exist", optname, arg)
        opts[[optname]] <- abspath(arg)
      }
  }
  
  # Prepare 2nd set of input functional filenames
  if (!is.null(opts$infuncs2)) {
      infiles2 <- sapply(as.character(read.table(opts$infuncs2)[,1]), function(fp) {
            if (!file.exists(fp))
                stop("One of the input functionals does not exist: ", fp)
            abspath(fp)
      })
      if (length(infiles2) != n) {
          vstop("Number of lines in %s doesn't match those in %s", 
                  opts$infuncs2, opts$infuncs1)
      }
  } else {
      infiles2 <- NULL
  }
  
  # Are we going to use a 2nd set of functionals?
  if (!is.null(infiles2) || !is.null(opts$brainmask2)) {
      use.set2 <- TRUE
      vcat(opts$verbose, "...using 2nd set of functional images")
      if (is.null(infiles2))
          infiles2 <- infiles1
  } else {
      use.set2 <- FALSE
      vcat(opts$verbose, "...NOT using 2nd set of functional images")
  }
  
  if (opts$debug) {
      verbosity <- 2
  } else if (opts$verbose) {
      verbosity <- 1
  } else {
      verbosity <- 0
  }
  
  if (opts$sparsedist) {
      method <- "icov"
  } else {
      method <- "pearson"
  }
  
  if (opts$sparseconn) {
      glasso <- TRUE
      if (!use.set2)
        stop("can only use --sparseconn when set --infiles2 and/or --brainmask2")
  } else {
      glasso <- FALSE
  }
  
  # design matrix
  if (!is.null(opts$regress)) {
      vcat(opts$verbose, "Reading in design matrix")
      tmp_fname <- opts$regress
      tmp <- as.matrix(read.table(opts$regress, header=TRUE))
      opts$regress <- big.matrix(nrow(tmp), ncol(tmp), type="double", shared=TRUE)
      opts$regress[,] <- tmp[,]; rm(tmp)
      k <- qlm_rank(opts$regress)
      if (k < ncol(opts$regress))
          vstop("design matrix (--regress %s) is rank deficient", tmp_fname)
      rm(tmp_fname)
  }
  
  
  ###
  # Read/Setup Masks
  ###
  
  vcat(opts$verbose, "Setting up masks")
  
  ## remove existing output
  if (opts$overwrite)
      stop("Right now the overwrite function isn't implemented")
  
  get_mask <- function(infiles, mask=NULL) {
      if (is.null(mask)) {
          hdr <- read.nifti.header(infiles[[1]])
          if (length(hdr$dim) == 2) {
              nvoxs <- hdr$dim[2]
          } else {
              nvoxs <- prod(hdr$dim[-length(hdr$dim)])
          }
          mask <- rep(TRUE, nvoxs)
      } else {
          mask <- read.mask(mask)
      }
      
      return(mask)
  }
  
  mask1 <- get_mask(infiles1, opts$brainmask1)
  if (use.set2) {
      mask2 <- get_mask(infiles2, opts$brainmask2)
  } else {
      mask2 <- NULL
  }
  
  
  ###
  # Set Memory Demands
  ###
  
  get_tpts <- function(x) {
      hdr <- read.nifti.header(x)
      n <- length(hdr$dim)
      if (n == 4) {
          return(hdr$dim[4])
      } else if (n == 2) {
          return(hdr$dim[1])
      } else {
          vstop("Input functional file '%s' must be 2 or 4 dimensions but is %i dimensional", x, n)
      }
  }
  
  nsubs <- length(infiles1)
  nvoxs1 <- sum(mask1)
  ntpts1 <- sapply(infiles1, get_tpts)
  if (use.set2) {
      nvoxs2 <- sum(mask2)
      ntpts2 <- sapply(infiles2, get_tpts)
      for (i in 1:nsubs) {
          if (ntpts1[i] != ntpts2[i]) {
              vstop("subject #%i does not have the same # of timepoints for the first and second functional datasets", i)
          }
      }
  } else {
      nvoxs2 <- NULL
  }
  opts <- get_subdist_memlimit(opts, nsubs, nvoxs1, ntpts1, nvoxs2)
  
  
  ###
  # Read/Prepare Functional Data
  ###
  
  vcat(opts$verbose, "Loading and masking functional data (Part 1)")
  ftype1 <- ifelse(opts$in2d1, "nifti2d", "nifti4d")
  reader1 <- gen_big_reader(ftype1, type="double", shared=parallel_forks)
  funclist1 <- load_and_mask_func_data2(infiles1, reader1, mask=mask1, 
                                        verbose=opts$verbose, scale=!glasso,  
                                        type="double", shared=parallel_forks)
  check1 <- check_func_data(infiles1[1], funclist1[1], extra=TRUE, 
                            verbose=opts$verbose, parallel=FALSE)
  check2 <- check_func_data(infiles1[-1], funclist1[-1], extra=opts$extrachecks, 
                            verbose=opts$verbose, parallel=parallel_forks)
  checks <- c(check1, check2)
  if (any(checks!=0)) {
      vcat(opts$verbose, "Bad data for following files:")
      vcat(opts$verbose, paste(infiles1[checks!=0], collapse="\n"))
      vstop("Quitting due to errors with 1st set of input functional data")
  }
  
  if (use.set2) {
      vcat(opts$verbose, "Loading and masking functional data (Part 2)")
      ftype2 <- ifelse(opts$in2d2, "nifti2d", "nifti4d")
      reader2 <- gen_big_reader(ftype2, type="double", shared=parallel_forks)
      funclist2 <- load_and_mask_func_data2(infiles2, reader2, mask=mask2, 
                                            verbose=opts$verbose, scale=!glasso,  
                                            type="double", shared=parallel_forks)
      check1 <- check_func_data(infiles2[1], funclist2[1], extra=opts$extrachecks, 
                                verbose=opts$verbose, parallel=parallel_forks)
      check2 <- check_func_data(infiles2[-1], funclist2[-1], extra=opts$extrachecks, 
                                verbose=opts$verbose, parallel=parallel_forks)
      checks <- c(check1, check2)
      if (any(checks!=0)) {
          vcat(opts$verbose, "Bad data for following files:")
          vcat(opts$verbose, paste(infiles2[checks!=0], collapse="\n"))
          vstop("Quitting due to errors with 2nd set of input functional data")
      }
  }
  
  invisible(gc(FALSE, TRUE))
  
  
  ###
  # Creating output directory
  ###
  
  vcat(opts$verbose, "Creating output directory and files")
  dists_list <- create_subdist(outdir, infiles1, mask1, infiles2, mask2, opts, shared=parallel_forks)  
  invisible(gc(FALSE, TRUE))
  
  end.time <- Sys.time()
  vcat(opts$verbose, "Setup is done. It took: %.2f minutes\n", 
       as.numeric(end.time-start.time, units="mins"))
  
  
  ###
  # Compute the subdists
  ###
  start.time <- Sys.time()
  
  vcat(opts$verbose, "Computing subject distances")
  if (use.set2) {
      checks <-  compute_subdist_wrapper3(funclist1, dists_list, 
                                        opts$blocksize, opts$superblocksize, 
                                        funclist2, 
                                        design_mat=opts$regress, 
                                        verbose=verbosity, parallel=parallel_forks, 
                                        ztransform=opts$ztransform, method=method, 
                                        glasso=glasso)
  } else {
      checks <- compute_subdist_wrapper(funclist1, dists_list, 
                                        opts$blocksize, opts$superblocksize, 
                                        design_mat=opts$regress, 
                                        verbose=verbosity, parallel=parallel_forks, 
                                        ztransform=opts$ztransform, method=method)
  }
  
  #vcat(opts$verbose, "...saving zchecks")  
  #hdr <- read.nifti.header(infiles1[1])
  #if (length(hdr$dim) == 4) {
  #    hdr$dim <- hdr$dim[1:3]; hdr$pixdim <- hdr$pixdim[1:3]
  #}
  #write.nifti(checks$sdist, hdr, mask1, odt="char", 
  #            outfile=file.path(outdir, "zcheck_subdist.nii.gz"))
  #write.nifti(checks$gdist, hdr, mask1, odt="char", 
  #            outfile=file.path(outdir, "zcheck_subdist_gower.nii.gz"))

  end.time <- Sys.time()
  vcat(opts$verbose, "Done! Total computation time: %.1f minutes\n", 
       as.numeric(end.time-start.time, units="mins"))  
  
}, warning = function(ex) {
  cat("\nA warning was detected: \n")
  cat(ex$message, "\n\n")
  cat("Called by: \n")
  print(ex$call)
  cat("\nSaving options...\n")
  save(saved_opts, file="called_options.rda")
}, error = function(ex) {
  cat("\nAn error was detected: \n")
  cat(ex$message, "\n\n")
  cat("Called by: \n")
  print(ex$call)
  cat("\nSaving options...\n")
  save(saved_opts, file="called_options.rda")
}, interrupt = function(ex) {
  cat("\nSaving options...\n")
  save(saved_opts, file="called_options.rda")
  cat("\nKill signal sent. Trying to clean up...\n")
  rm(list=ls())
  gc(FALSE)
  cat("...success\n")
}, finally = {
  cat("\nRemoving everything from memory\n")
  rm(list=ls())
  gc(FALSE)
  cat("...sucesss\n")
})
