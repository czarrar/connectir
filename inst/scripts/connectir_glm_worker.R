suppressPackageStartupMessages(library("optparse"))

# Make option list
option_list <- list(
    make_option(c("-i", "--infuncs1"), type="character", default=NULL, dest="infuncs1", help="File containing paths of different functional images (nifti or text files) in one column. Each node/voxel in these images will act as a 'seed'. (required)", metavar="file"),
    make_option("--brainmask1", type="character", default=NULL, help="When computing each whole-brain connectivity map, this mask will restrict which parts of the whole-brain are to be considered", metavar="file", dest="brainmask1"), 
    make_option(c("-r", "--regressors"), type="character", default=NULL, dest="evs", help="File containing regressors (required)", metavar="file"), 
    make_option(c("-n", "--contrasts"), type="character", default=NULL, dest="contrasts", help="File containing contrasts (required)", metavar="file"), 
    make_option("--infuncs2", type="character", default=NULL, dest="infuncs2", help="File containing paths of different functional images (nifti or text files) in one column. Each node/voxel in these images will act as a 'target', that is voxels from -i/--infuncs will be correlated with voxels in these images. The default is to have these images be the same as those specified in -i/--infuncs. (optional)", metavar="file"),
    make_option("--brainmask2", type="character", default=NULL, help="When computing each whole-brain connectivity map, this mask will restrict which parts of the whole-brain are to be considered only for --infuncs2", metavar="file", dest="brainmask2"), 
    make_option("--ztransform", action="store_true", default=FALSE, dest="ztransform", help="Fischer Z-Transform the correlations before relating the regressors to connectivity"), 
    make_option("--summarize", action="store_true", default=FALSE, dest="summarize", help="Output summary measures of t-tests"), 
    make_option("--memlimit", type="double", default=4, dest="memlimit", help="If blocksize is set to auto (--blocksize=0), then will set the blocksize to use a maximum of RAM specified by this option  [default: %default]", metavar="number"),
    make_option(c("-c", "--forks"), type="integer", default=1, help="Number of computer processors to use in parallel by forking the complete processing stream [default: %default]", metavar="number"),
    make_option(c("-t", "--threads"), type="integer", default=1, help="Number of computer processors to use in parallel by multi-threading matrix algebra operations [default: %default]", metavar="number"),
    make_option("--extrachecks", action="store_true", default=FALSE, help="Will do a more rigorous check of the input functionals before any calculations"),
    make_option("--overwrite", action="store_true", default=FALSE, help="Overwrite output that already exists (default is not to overwrite already existing output)"), 
    make_option(c("-q", "--quiet"), action="store_false", dest="verbose", help="Print little output"), 
    make_option(c("-v", "--verbose"), action="store_true", default=TRUE, help="Print extra output [default]")
)

# Make class/usage
parser <- OptionParser(usage = "%prog [options] output-directory", 
                       option_list=option_list, add_help_option=TRUE)

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
  parallel <- ifelse(opts$forks == 1, FALSE, TRUE)

  ###
  # Check/Setup Required Inputs
  ###
  
  # Check required
  vcat(opts$verbose, "Checking user inputs")
  outdir <- abspath(args[1])
  if (file.exists(outdir) && !opts$overwrite)
      stop("Output '", outdir, "' already exists, you can use --overwrite")
  if (is.null(opts$infuncs1))
      stop("You must specify the -i/--infuncs option")
  if (is.null(opts$evs))
      stop("You must specify the -e/--evs option")
  if (is.null(opts$contrasts))
      stop("You must specify the -c/--contrasts option")
  if (!file.exists(opts$infuncs1))
      stop("The file specified by -i/--infuncs does not exist")
  if (!file.exists(opts$evs))
      vstop("The regressors file '%s' does not exist", opts$evs)
  if (!file.exists(opts$contrasts))
      vstop("The contrast file '%s' does not exist", opts$contrasts)
  if (!is.null(opts$infuncs2) && !file.exists(opts$infuncs2))
      stop("The file specified by --infuncs2 does not exist")
  
  func_files1 <- as.character(read.table(opts$infuncs1)[,1])
  if (!is.null(opts$infuncs2)) {
      func_files2 <- as.character(read.table(opts$infuncs2)[,1])
  } else {
      func_files2 <- NULL
  }
  
  # Running glm wrapper
  wrap_glm(func_files1, opts$brainmask1, opts$evs, opts$contrasts, 
           outdir, summarize=opts$summarize, overwrite=opts$overwrite, 
           func_files2=func_files2, mask_file2=opts$brainmask2, 
           memlimit=opts$memlimit, extra_checks=opts$extrachecks, 
           verbose=opts$verbose, parallel=parallel, shared=parallel, 
           ztransform=opts$ztransform) 

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
