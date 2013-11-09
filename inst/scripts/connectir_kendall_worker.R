suppressPackageStartupMessages(library("optparse"))

# Make option list
option_list <- list(
    make_option(c("-i", "--infuncs1"), type="character", default=NULL, dest="infuncs1", help="File containing paths of different functional images (nifti or text files) in one column. Each node/voxel in these images will act as a 'seed'. (required)", metavar="file"),
    make_option("--infuncs2", type="character", default=NULL, dest="infuncs2", help="File containing paths of different functional images (nifti or text files) in one column. Each node/voxel in these images will act as a 'target', that is voxels from -i/--infuncs will be correlated with voxels in these images. The default is to have these images be the same as those specified in -i/--infuncs. (optional)", metavar="file"),
    make_option("--ztransform", action="store_true", default=FALSE, dest="ztransform", help="Fischer Z-Transform the correlations before calculating the similarity of connectivity maps across participants"),
    make_option("--brainmask1", type="character", default=NULL, help="When computing each whole-brain connectivity map, this mask will restrict which parts of the whole-brain are to be considered", metavar="file", dest="brainmask1"),
    make_option("--brainmask2", type="character", default=NULL, help="When computing each whole-brain connectivity map, this mask will restrict which parts of the whole-brain are to be considered only for --infuncs2", metavar="file", dest="brainmask2"), 
    make_option("--regress", type="character", default=NULL, help="A design matrix (space delimeted file where first row is a header) containing variables to regress out of each voxel's whole-brain connectivity maps before comparing distances between subjects", metavar="file"), 
    make_option("--memlimit", type="double", default=4, dest="memlimit", help="If blocksize is set to auto (--blocksize=0), then will set the blocksize to use a maximum of RAM specified by this option  [default: %default]", metavar="number"),
    make_option(c("-c", "--forks"), type="integer", default=1, help="Number of computer processors to use in parallel by forking the complete processing stream [default: %default]", metavar="number"),
    make_option(c("-t", "--threads"), type="integer", default=1, help="Number of computer processors to use in parallel by multi-threading matrix algebra operations [default: %default]", metavar="number"),
    make_option("--extrachecks", action="store_true", default=FALSE, help="Will do a more rigorous check of the input functionals before any calculations"),
    make_option("--overwrite", action="store_true", default=FALSE, help="Overwrite output that already exists (default is not to overwrite already existing output)"), 
    make_option(c("-q", "--quiet"), action="store_false", dest="verbose", help="Print little output"), 
    make_option(c("-v", "--verbose"), action="store_true", default=TRUE, help="Print extra output [default]")
)

# Make class/usage
parser <- OptionParser(usage = "%prog [options] output", option_list=option_list, add_help_option=TRUE)

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
  output <- abspath(args[1])
  if (file.exists(output) && !opts$overwrite)
      stop("Output '", output, "' already exists, you can use --overwrite")
  if (is.null(opts$infuncs1))
      stop("You must specify the -i/--infuncs option")
  if (!file.exists(opts$infuncs1))
      stop("The file specified by -i/--infuncs does not exist")
  if (!is.null(opts$infuncs2) && !file.exists(opts$infuncs2))
      stop("The file specified by --infuncs2 does not exist")
  if (!is.null(opts$regress) && !file.exists(opts$regress))
      stop("The file specified by --regress does not exist")
  
  func_files1 <- as.character(read.table(opts$infuncs1)[,1])
  if (!is.null(opts$infuncs2)) {
      func_files2 <- as.character(read.table(opts$infuncs2)[,1])
  } else {
      func_files2 <- NULL
  }
    
  # Running kendall wrapper
  wrap_kendall(func_files1, opts$brainmask1, 
               func_files2, opts$brainmask2, 
               design_mat=opts$regress, 
               out_file=output, overwrite=opts$overwrite, 
               verbose=opts$verbose, parallel=parallel, shared=parallel, 
               memlimit=opts$memlimit, extra_checks=opts$extrachecks, 
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
