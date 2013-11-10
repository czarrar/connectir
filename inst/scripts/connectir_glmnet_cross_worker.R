suppressPackageStartupMessages(library("optparse"))

# Make option list
option_list <- list(
    make_option(c("-i", "--subdist"), type="character", default=NULL, help="Input subject distances descriptor (*.desc) file (required)", metavar="file"), 
    make_option(c("-m", "--mask"), type="character", default=NULL, help="Brain mask file (required)", metavar="file"), 
    make_option(c("-l", "--labels"), type="character", default=NULL, help="File with labels/responses where # of rows correspond to # of subjects in the subject distances matrices (required)", metavar="file"), 
    make_option(c("-c", "--forks"), type="integer", default=1, help="Number of computer processors to use in parallel by forking the complete processing stream [default: %default]", metavar="number"),
    make_option(c("-t", "--threads"), type="integer", default=1, help="Number of computer processors to use in parallel by multi-threading matrix algebra operations [default: %default]", metavar="number"),
    make_option("--cross", type="integer", default=10, help="Number of folds for cross-validation (default: %default)", metavar="option"),  
    make_option("--family", type="character", default=NULL, help="Type of classification: 'gaussian', 'binomial', 'poisson', 'multinomial', 'cox' (default is to auto-select between binomial, multinomial, and gaussian)", metavar="option"), 
    make_option("--memlimit", type="double", default=1, dest="memlimit", help="Maximum amount of RAM to use. It is preferable to keep this number as small as possible for speed reasons (this rule of thumb just applies to this script and connectir_kmeans_cross).  [default: %default]", metavar="RAM"),
    make_option("--overwrite", action="store_true", default=FALSE, help="Overwrite output if it already exists (default is not to overwrite already existing output)"),
    make_option(c("-v", "--verbose"), action="store_true", default=TRUE, help="Print extra output [default]"),
    make_option(c("-q", "--quiet"), action="store_false", dest="verbose", help="Print little output")
)

# Make class/usage
parser <- OptionParser(usage = "%prog [options] outfile", 
                       option_list=option_list, add_help_option=TRUE)

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
  if (is.null(opts$subdist))
      stop("Must specify -i/--subdist")
  if (is.null(opts$mask))
      stop("Must specify -m/--mask")
  if (is.null(opts$labels))
      stop("Must specify -l/--labels")
  if (getext(opts$subdist) != "desc")
      stop("Subject distances file (-i/--subdist) must have a '.desc' extension")
  opts$prefix <- args[1]
  
  
  ###
  # Compute Baby Compute
  ###
  start.time <- Sys.time()
  wrap_glmnet_subdist_cross(opts$subdist, opts$mask, opts$labels, 
                            out_prefix=opts$prefix, overwrite=opts$overwrite, 
                            family=opts$family, standardize=TRUE, cross=opts$cross, 
                            memlimit=opts$memlimit, parallel=parallel_forks, 
                            verbose=opts$verbose)
  end.time <- Sys.time()
  vcat(opts$verbose, "SVM is done! It took: %.2f minutes\n", 
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

  
