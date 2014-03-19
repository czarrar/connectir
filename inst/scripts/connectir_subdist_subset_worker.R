suppressPackageStartupMessages(library("optparse"))


# Make option list
option_list <- list(
    make_option(c("-i", "--input"), type="character", default=NULL, help="Subject distances (use the *.desc extension) (required)", metavar="file"),
    make_option(c("-v", "--voxels"), type="character", default=NULL, help="Range of voxels (columns) to extract. For example, '1' or '1:10' so basically it can be any R expression. Need this or -m/--mask", metavar="number"),
    make_option(c("-m", "--mask"), type="character", default=NULL, help="Mask indicating voxels (columns) to extract with value of 2 and other voxels used for subject distances with value 1. Need this or -v/--voxels", metavar="file"),
    make_option(c("-c", "--forks"), type="integer", default=1, help="Number of computer processors to use in parallel by forking the complete processing stream [default: %default]", metavar="number"),
    make_option(c("-t", "--threads"), type="integer", default=1, help="Number of computer processors to use in parallel by multi-threading matrix algebra operations [default: %default]", metavar="number"),
    make_option("--blocksize", type="integer", default=0, dest="blocksize", help="How many sets of columns or voxels should used in each iteration of computing the element-wise average (0 = auto) [default: %default]", metavar="number"),
    make_option("--memlimit", type="double", default=4, dest="memlimit", help="If blocksize is set to auto (--blocksize=0), then will set the blocksize to use a maximum of RAM specified by this option  [default: %default]", metavar="number"),
    make_option("--overwrite", action="store_true", default=FALSE, help="Overwrite output if it already exists (default is not to overwrite already existing output)"),
    make_option("--debug", action="store_true", default=FALSE, help="Like verbose but will also print more helpful error messages when --forks is >1"), 
    make_option("--verbose", action="store_true", default=TRUE, help="Print extra output [default]"),
    make_option("--quiet", action="store_false", dest="verbose", help="Print little output")
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
  if (is.null(opts$input))
      stop("must specify -i/--input")
  if (opts$overwrite)
      stop("Right now the overwrite function isn't implemented")  
  if (is.null(opts$voxels) && is.null(opts$mask))
      stop("Must specify either -v/--voxels or -m/--mask")
  if (!is.null(opts$voxels) && !is.null(opts$mask))
      stop("Must specify either -v/--voxels or -m/--mask (not both)")
  
  ###
  # Setup inputs
  ###
  vcat(opts$verbose, "Setting/checking inputs")
  
  # check/set input
  opts$input <- abspath(opts$input)
  if (!file.exists(opts$input))
      vstop("-i/--input '%s' does not exist", opts$input)
  
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
  
  vcat(opts$verbose, "Reading in input")
  
  imat <- attach.big.matrix(opts$input)
  
  if (!is.null(opts$mask)) {
      mask <- read.nifti.image(opts$mask)
      # check that mask only has needed values
      if (!all(sort(unique(mask)) != c(0,1,2))) {
          vstop("Mask should only have 0,1,2 as it's unique values. Something else was found")
      }
      # check that mask is same size as imat
      if (sum(mask>0) != ncol(imat)) {
          vstop("mask should have the same number of non-zero values as columns in the subject distances")
      }
      voxs <- which(mask==2)
   } else if (!is.null(opts$voxels)) {
       voxs <- eval(parse(text=opts$voxels))
   }
  
  
  ###
  # Setup outputs
  ###
  
  vcat(opts$verbose, "Setup outputs")
  
  omat <- big.matrix(nrow(imat), length(voxs), type=typeof(imat), init=0, 
                     backingfile=paste(outprefix, ".bin", sep=""), 
                     descriptorfile=paste(outprefix, ".desc", sep=""),
                     backingpath=dirname(outprefix))
  
  ###
  # Subset / Save
  ###
  vcat(opts$verbose, "Getting subset of voxels")
  
  deepcopy(x=imat, cols=voxs, y=omat)
        
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
