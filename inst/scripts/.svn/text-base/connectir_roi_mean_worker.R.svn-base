suppressPackageStartupMessages(library("optparse"))

# Make option list
option_list <- list(
    make_option(c("-i", "--input"), type="character", default=NULL, help="Input 4D functional image (required)", metavar="file"),
    make_option(c("-r", "--roi"), type="character", default=NULL, help="3D image with ROIs (required)", metavar="file"),
    make_option(c("-m", "--mask"), type="character", default=NULL, help="3D brain mask (optional)", metavar="file"),
    make_option("--outtype", type="character", default="nifti", help="Type of output: nifti or text (default: %default)", metavar="option"),
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
  
  # set output
  opts$outfile <- args[1]
  if (opts$outtype == "text") {
      opts$outfile <- paste(rmext(opts$outfile), "txt", sep=".")
  } else if (opts$outtype == "nifti") {
      opts$outfile <- paste(rmext(opts$outfile), "nii.gz", sep=".")
  } else {
      vstop("unrecognized output type '%s'; must be text or nifti", opts$outtype)
  }
  
  
  ###
  # Check Inputs
  ###
  vcat(opts$verbose, "Checking options")
  
  if (is.null(opts$input))
      stop("Must specify input (-i/--input)")
  if (is.null(opts$roi))
      stop("Must specify roi (-r/--roi)")
  
  if (!file.exists(opts$input))
      vstop("Input '%s' doesn't exist", opts$input)
  if (!file.exists(opts$roi))
      vstop("ROI file '%s' doesn't exist", opts$roi)
  if (!is.null(opts$mask) && !file.exists(opts$mask))
      vstop("Mask '%s' doesn't exist", opts$mask)
  if (file.exists(opts$outfile) && !opts$overwrite)
      vstop("Output '%s' already exists (consider using --overwrite)", opts$outfile)
  
  
  ###
  # Setup
  ###
  vcat(opts$verbose, "Setup inputs")
  
  input <- abspath(opts$input)
  roi <- abspath(opts$roi)
  if (!is.null(opts$mask))
      mask <- abspath(opts$mask)
  else
      mask <- NULL
  outfile <- abspath(opts$outfile)
  outtype <- opts$outtype
  overwrite <- opts$overwrite
  verbose <- opts$verbose
  
  
  ###
  # Run IT!
  ###
  
  start.time <- Sys.time()
  roi_mean_wrapper(input, roi, mask, 
                   out_file=outfile, outtype=outtype, overwrite=overwrite, 
                   verbose=verbose)
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
