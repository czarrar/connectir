suppressPackageStartupMessages(library("optparse"))


# Make option list
option_list <- list(
    make_option(c("-a", "--amat"), type="character", default=NULL, help="First file-backed big matrix such as the subject distances (use the *.desc extension) (required)", metavar="file"),
    make_option(c("-b", "--bmat"), type="character", default=NULL, help="Second file-backed big matrix (use the *desc extension) (required)", metavar="file"),
    make_option(c("-c", "--forks"), type="integer", default=1, help="Number of computer processors to use in parallel by forking the complete processing stream [default: %default]", metavar="number"),
    make_option(c("-t", "--threads"), type="integer", default=1, help="Number of computer processors to use in parallel by multi-threading matrix algebra operations [default: %default]", metavar="number"),
    make_option("--blocksize", type="integer", default=0, dest="blocksize", help="How many sets of columns or voxels should used in each iteration of computing the element-wise average (0 = auto) [default: %default]", metavar="number"),
    make_option("--memlimit", type="double", default=4, dest="memlimit", help="If blocksize is set to auto (--blocksize=0), then will set the blocksize to use a maximum of RAM specified by this option  [default: %default]", metavar="number"),
    make_option("--overwrite", action="store_true", default=FALSE, help="Overwrite output if it already exists (default is not to overwrite already existing output)"),
    make_option(c("-d", "--debug"), action="store_true", default=FALSE, help="Like verbose but will also print more helpful error messages when --forks is >1"), 
    make_option(c("-v", "--verbose"), action="store_true", default=TRUE, help="Print extra output [default]"),
    make_option(c("-q", "--quiet"), action="store_false", dest="verbose", help="Print little output")
)

# Make class/usage
parser <- OptionParser(usage = "%prog [options] outprefix (extension should not exist for outprefix)", 
                       option_list=option_list, 
                       add_help_option=TRUE)

# Parse
parser_out <- parse_args(parser, positional_arguments = TRUE)
args <- parser_out$args
opts <- parser_out$options

# Check options/arguments
if (length(args) != 1) {
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
  # Check Inputs
  ###
  vcat(opts$verbose, "Checking options")
  
  # check variables
  if (is.null(opts$amat) || is.null(opts$bmat))
      stop("must specify --amat or --bmat")
  if (opts$overwrite)
      stop("Right now the overwrite function isn't implemented")  
  
  ###
  # Setup inputs
  ###
  vcat(opts$verbose, "Setting/checking inputs")
  
  # check/set input
  opts$amat <- abspath(opts$amat)
  opts$bmat <- abspath(opts$bmat)
  if (!file.exists(opts$amat))
      vstop("-a/--amat '%s' does not exist", opts$amat)
  if (!file.exists(opts$bmat))
      vstop("-b/--bmat '%s' does not exist", opts$bmat)
  
  # check/set output
  outprefix <- abspath(args[1])
  if (!file.exists(dirname(args[1]))) {
      vcat(opts$verbose, "Output directory '%s' does not exist, creating", 
           dirname(args[1]))
      dir.create(dirname(args[1]))
  }
  
  if (opts$debug) {
      verbosity <- 2
  } else if (opts$verbose) {
      verbosity <- 1
  } else {
      verbosity <- 0
  }
  
  ###
  # Read inputs
  ###
  
  vcat(opts$verbose, "Reading in inputs")
  
  amat <- attach.big.matrix(opts$amat)
  bmat <- attach.big.matrix(opts$bmat)
  
  # check dimensions
  nr <- nrow(amat); nc <- ncol(amat); type <- typeof(amat)
  if ((nr != nrow(bmat)) || (nc != ncol(bmat))) {
      vstop("Dimensions for -a/--amat (%i,%i) and -b/--bmat (%i,%i) do not match", nr, nc, nrow(bmat), ncol(bmat))
  }
  if (type != typeof(bmat)) {
      vstop("-a/--amat (%s) and -b/--bmat (%s) do not have the same data types", type, typeof(bmat))
  }
  
  ###
  # Setup outputs
  ###
  
  vcat(opts$verbose, "Setup outputs")
  
  omat <- big.matrix(nr, nc, type=type, init=0, 
                     backingfile=paste(outprefix, ".bin", sep=""), 
                     descriptorfile=paste(outprefix, ".desc", sep=""),
                     backingpath=dirname(outprefix))
  
  ###
  # Average / Save
  ###
  vcat(opts$verbose, "Averaging")
  
  daxpy(ALPHA=0.5, X=amat, Y=omat)
  daxpy(ALPHA=0.5, X=bmat, Y=omat)
      
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
