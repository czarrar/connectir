suppressPackageStartupMessages(library("optparse"))


# Make option list
option_list <- list(
    make_option(c("-i", "--subdist"), type="character", default=NULL, help="Input subject distances descriptor file (required)", metavar="file"),
    make_option(c("-m", "--model"), type="character", default=NULL, help="Filename of a comma separated file with participant info in R friendly format where column names correspond to formula values... (required)", metavar="csv"),
    make_option("--whichsubs", type="character", default=NULL, help="Filename with a list of subject indices to use from the subject distance matrices (default is to use all of them)", metavar="text-file"),
    make_option("--expr", type="character", default=NULL, help="An expression based on the model that is used to restrict the subjects examined (can either use this or --whichsubs, not both)", metavar="expression"),
    make_option("--nogower", action="store_true", default=FALSE, help="Don't run and save gower centered distance matrices"), 
    make_option(c("-c", "--forks"), type="integer", default=1, help="Number of computer processors to use in parallel by forking the complete processing stream [default: %default]", metavar="number"),
    make_option(c("-t", "--threads"), type="integer", default=1, help="Number of computer processors to use in parallel by multi-threading matrix algebra operations [default: %default]", metavar="number"),
    make_option("--blocksize", type="integer", default=0, dest="blocksize", help="How many sets of voxels should used in each iteration of computing the pseudo F-statistics (0 = auto) [default: %default]", metavar="number"),
    make_option("--memlimit", type="double", default=4, dest="memlimit", help="If blocksize is set to auto (--blocksize=0), then will set the blocksize to use a maximum of RAM specified by this option  [default: %default]", metavar="number"),
    make_option("--overwrite", action="store_true", default=FALSE, help="Overwrite output if it already exists (default is not to overwrite already existing output)"),
    make_option(c("-d", "--debug"), action="store_true", default=FALSE, help="Like verbose but will also print more helpful error messages when --forks is >1"), 
    make_option(c("-v", "--verbose"), action="store_true", default=TRUE, help="Print extra output [default]"),
    make_option(c("-q", "--quiet"), action="store_false", dest="verbose", help="Print little output")
)
    
# Make class/usage
parser <- OptionParser(usage = "%prog [options] outprefix", 
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
  if (is.null(opts$subdist))
      stop("must specify -i/--subdist option")
  if (is.null(opts$model))
      stop("must specify -m/--model option")
  if (opts$overwrite)
      stop("Right now the overwrite function isn't implemented")
  if (!is.null(opts$expr) && !is.null(opts$whichsubs))
      stop("cannot specify both --expr and --whichsubs")
  
  
  ###
  # Setup inputs
  ###
  vcat(opts$verbose, "Setting/checking inputs")
  
  # check/set input
  opts$subdist <- abspath(opts$subdist)
  if (!file.exists(opts$subdist))
      vstop("-i/--subdist '%s' does not exist", opts$subdist)
  if (!file.exists(opts$model))
      vstop("-m/--model '%s' does not exist", opts$model)
  
  # check/set output
  if (!file.exists(dirname(args[1]))) {
      vcat(opts$verbose, "Output directory '%s' does not exist, creating", 
           dirname(args[1]))
      dir.create(dirname(args[1]))
  }
  outprefix <- abspath(args[1])
  ## subject distances
  s.path <- list()
  s.path$bpath <- dirname(outprefix)
  s.path$bfile <- sprintf("%s_subdist.bin", basename(outprefix))
  s.path$dfile <- sprintf("%s_subdist.desc", basename(outprefix))
  outfile <- file.path(s.path$bpath, s.path$dfile)
  if (file.exists(outfile))
      vstop("Output subject distances '%s' already exists", outfile)
  ## gower distances
  g.path <- list()
  g.path$bpath <- dirname(outprefix)
  g.path$bfile <- sprintf("%s_subdist_gower.bin", basename(outprefix))
  g.path$dfile <- sprintf("%s_subdist_gower.desc", basename(outprefix))
  outfile <- file.path(g.path$bpath, g.path$dfile)
  if (!opts$nogower && file.exists(outfile))
      vstop("Output gower distances '%s' already exists", outfile)
  ## model
  m.path <- sprintf("%s_model.csv", outprefix)
  
  if (opts$debug) {
      verbosity <- 2
  } else if (opts$verbose) {
      verbosity <- 1
  } else {
      verbosity <- 0
  }
  
  # data dimensions
  vcat(opts$verbose, "...data dimensions")
  tmp <- attach.big.matrix(opts$subdist)
  nsubs <- sqrt(nrow(tmp))
  nvoxs <- ncol(tmp)
  rm(tmp); gc(FALSE, TRUE)
  
  # model
  vcat(opts$verbose, "...model")
  model <- read.csv(opts$model)
  
  # filter subjects
  if (!is.null(opts$whichsubs)) {
      vcat(opts$verbose, "...whichsubs")
      which.subs <- as.numeric(read.table(opts$whichsubs)[,1])
      if (all(which.subs==0 | which.subs==1))
          which.subs <- which(which.subs==1)
      if (length(which.subs) == 0)
          stop("no subjects left to analyze based on --whichsubs")
  } else if (!is.null(opts$expr)) {
      vcat(opts$verbose, "...expr")
      which.subs <- eval(parse(text=sprintf("with(model, which(%s))", opts$expr)))
      if (length(which.subs) == 0)
          stop("no subjects left to analyze based on --expr")
  } else {
      which.subs <- 1:nrow(model)
  }
  ## model
  if (nrow(model) == nsubs) {
      model <- model[which.subs,]
  } else if (nrow(model) != length(which.subs)) {
      stop(paste("# of rows in model don't match # of subjects in", 
                 "distance matrix or # of filtered subjects"))
  } else {
      stop("error")
  }
  
  ###
  # Filter / Save
  ###
  vcat(opts$verbose, "Filtering and/or Gowering")
  
  # convert
  xdist <- filter_subdist_fb(opts$subdist, which.subs, s.path, opts$memlimit, 
                             parallel=parallel_forks, verbose=opts$verbose, 
                             gower=!opts$nogower, gower.paths=g.path)
  
  # check
  nsubs <- sqrt(nrow(xdist))
  vcat(opts$verbose, "...checking input")
  tmp <- matrix(xdist[,1], nsubs, nsubs)
  check_dmat(tmp)
  rm(tmp); invisible(gc(FALSE, TRUE))
  
  vcat(opts$verbose, "Saving model")
  write.csv(model, file=m.path)
    
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

  
  