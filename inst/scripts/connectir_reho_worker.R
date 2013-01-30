suppressPackageStartupMessages(library("optparse"))


# Make option list
option_list <- list(
    make_option(c("-i", "--input"), type="character", default=NULL, help="Input 4D functional image (required)", metavar="file"),
    make_option(c("-m", "--mask"), type="character", default=NULL, help="3D brain mask (required)", metavar="file"),
    make_option("--min", type="double", default=0.5, help="Percent of neighboring voxels that must be non-zero if voxel is to be used (otherwise will have a value o 0) [default: %default]", metavar="0-1"),
    make_option("--nei1", type="integer", default=1, help="voxel distance #1 (more of an internal option) [default: %default]", metavar="positive-number"),
    make_option("--nei2", type="integer", default=3, help="voxel distance #2 (more of an internal option) [default: %default]", metavar="positive-number"),
    make_option(c("-c", "--forks"), type="integer", default=1, help="Number of computer processors to use in parallel by forking the complete processing stream [default: %default]", metavar="number"),
    make_option("--overwrite", action="store_true", default=FALSE, help="Overwrite output if it already exists"),
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
  set_parallel_procs(opts$forks, 1, opts$verbose)  
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
  
  overwrite <- opts$overwrite
  verbose <- opts$verbose
  nei <- opts$nei1
  nei.dist <- opts$nei2
  min.nei <- opts$min
  
  
  ###
  # Run IT!
  ###
  
  start.time <- Sys.time()
                     
  rs <- wrap_reho(input, mask, outfile, 
                  overwrite=overwrite, verbose=verbose, parallel=parallel, 
                  nei=nei, nei.dist=nei.dist, min.nei=min.nei)  
  
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
