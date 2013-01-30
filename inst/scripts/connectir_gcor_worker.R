suppressPackageStartupMessages(library("optparse"))


# Make option list
option_list <- list(
    make_option(c("-i", "--input"), type="character", default=NULL, help="Input 4D functional image (required)", metavar="file"),
    make_option(c("-m", "--mask"), type="character", default=NULL, help="3D brain mask (required)", metavar="file"),
    make_option("--ztransform", action="store_true", default=FALSE, help="Z-transform the data => average => reverse transform"),
    make_option("--thresh", type="double", default=NULL, help="Threshold to set (note: --type must also be set for this to work)", metavar="number"),
    make_option("--type", type="character", default=NULL, help="Type of threshold: positive (r>i) or negative (r<i)", metavar="option"),
    make_option("--bin", action="store_true", default=FALSE, help="Binarize the correlations based on a threshold and then take the sum; in other words, how many connections are greater than or less than a specific threshold (must set --type as well)"),
    make_option(c("-c", "--forks"), type="integer", default=1, help="Number of computer processors to use in parallel by forking the complete processing stream [default: %default]", metavar="number"),
    make_option(c("-t", "--threads"), type="integer", default=1, help="Number of computer processors to use in parallel by multi-threading matrix algebra operations [default: %default]", metavar="number"),
    make_option("--blocksize", type="integer", default=0, dest="blocksize", help="How many sets of voxels should used in each iteration of computing the pseudo F-statistics (0 = auto) [default: %default]", metavar="number"),
    make_option("--memlimit", type="double", default=4, dest="memlimit", help="If blocksize is set to auto (--blocksize=0), then will set the blocksize to use a maximum of RAM specified by this option  [default: %default]", metavar="number"),
    make_option("--overwrite", action="store_true", default=FALSE, help="Overwrite output if it already exists"),
    make_option(c("-d", "--debug"), action="store_true", default=FALSE, help="Like verbose but will also print more helpful error messages when --forks is >1"), 
    make_option(c("-v", "--verbose"), action="store_true", default=TRUE, help="Print extra output [default]"),
    make_option(c("-q", "--quiet"), action="store_false", dest="verbose", help="Print little output")
)

# Make class/usage
parser <- OptionParser(usage = "%prog [options] outfile", 
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
  
  # set output
  opts$outfile <- args[1]
  
  
  ###
  # Check Inputs
  ###
  vcat(opts$verbose, "Checking options")
  
  if (is.null(opts$input))
      stop("Must specify input (-i/--input)")
  if (is.null(opts$mask))
      stop("Must specify mask (-m/--mask)")
  if (is.null(opts$thresh) && !is.null(opts$type))
      stop("If setting thresh, must specify the type of threshold")
  if (!is.null(opts$thresh) && is.null(opts$type))
      stop("If setting type of threshold, must specify threshold")  
  if (opts$bin == T && (is.null(opts$type) || is.null(opts$thresh)))
      stop("If using --bin, must set --type and --thresh")
  
  if (!file.exists(opts$input))
      vstop("Input '%s' doesn't exist", opts$input)
  if (!file.exists(opts$mask))
      vstop("Mask '%s' doesn't exist", opts$mask)
  if (file.exists(opts$outfile) && !opts$overwrite)
      vstop("Output '%s' already exists (consider using --overwrite)", opts$outfile)
  
  ###
  # Setup
  ###
  vcat(opts$verbose, "Setup inputs")
  
  parallel <- opts$forks > 1
  
  input <- abspath(opts$input)
  mask <- abspath(opts$mask)
  outfile <- abspath(opts$outfile)
  
  if (is.null(opts$thresh))
      thresh <- 0
  else
      thresh <- opts$thresh
  
  if (is.null(opts$type)) {
      threshType <- 0
  } else if (opts$type == "positive") {
      threshType <- 1
  } else if (opts$type == "negative") {
      threshType <- 3
  } else {
      vstop("Unrecognized --type '%s'", opts$type)
  }
  
  if (opts$bin)
      threshType <- threshType + 1
  
  ztransform <- opts$ztransform
  overwrite <- opts$overwrite
  memlimit <- opts$memlimit
  blocksize <- opts$blocksize
  verbose <- opts$verbose
  
  ###
  # Run IT!
  ###
  
  start.time <- Sys.time()
  
  gs <- wrap_gcor(input, mask, outfile, 
                  blocksize=blocksize, memlimit=memlimit, 
                  overwrite=overwrite, verbose=verbose, parallel=parallel, 
                  ztransform=ztransform, thresh=thresh, threshType=threshType)
  
  end.time <- Sys.time()
  vcat(verbose, "Done! Total computation time: %.1f minutes\n", 
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
