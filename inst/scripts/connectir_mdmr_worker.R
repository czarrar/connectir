suppressPackageStartupMessages(library("optparse"))


# Make option list
option_list <- list(
    make_option(c("-i", "--indir"), type="character", default=NULL, help="Input subdist directory (required)", metavar="directory"),
    make_option("--subdist", type="character", default=NULL, help="Input subject distances descriptor file (optional, defaults to one in --indir)", metavar="file"),
    make_option(c("-f", "--formula"), type="character", default=NULL, help="a typical R model formula that specifies the factors or continuous variables that may expain the variance in each voxel's subject distance matrix", metavar="'A + B:C'"),
    make_option(c("-m", "--model"), type="character", default=NULL, help="Filename of a comma separated file with participant info in R friendly format where column names correspond to formula values... (required)", metavar="csv"),
    make_option("--strata", type="character", default=NULL, help="Only compute permutations within groups, you can specify the name of a column in your '--model' that indicates these groups (optional)", metavar="name"),
    make_option(c("-p", "--permutations"), type="integer", default=4999, help="Number of permutations to conduct for each voxel [default: %default]", metavar="number"),
    make_option("--factors2perm", type="character", default=NULL, help="Which factors (e.g., A and B) to permute from the formula specified [default: all of them]", metavar="'A,B'"),
    make_option(c("-c", "--forks"), type="integer", default=1, help="Number of computer processors to use in parallel by forking the complete processing stream [default: %default]", metavar="number"),
    make_option(c("-t", "--threads"), type="integer", default=1, help="Number of computer processors to use in parallel by multi-threading matrix algebra operations [default: %default]", metavar="number"),
    make_option(c("-j", "--jobs"), type="integer", default=NULL, help="Number of SGE jobs to submit (if not set, won't use SGE)", metavar="number"), 
    make_option("--blocksize", type="integer", default=0, dest="blocksize", help="How many sets of voxels should used in eaciteration of computing the pseudo F-statistics (0 = auto) [default: %default]", metavar="number"),
    make_option("--memlimit", type="double", default=6, dest="memlimit", help="If blocksize is set to auto (--blocksize=0), then will set the blocksize to use a maximum of RAM specified by this option  [default: %default]", metavar="number"),
    make_option("--save-perms", action="store_true", default=FALSE, dest="saveperms", help="Save all the permuted psuedo-F stats? [default: %default]"),
    make_option("--overwrite", action="store_true", default=FALSE, help="Overwrite output if it already exists (default is not to overwrite already existing output)"),
    make_option(c("-d", "--debug"), action="store_true", default=FALSE, help="Like verbose but will also print more helpful error messages when --forks is >1"), 
    make_option(c("-v", "--verbose"), action="store_true", default=TRUE, help="Print extra output [default]"),
    make_option(c("-q", "--quiet"), action="store_false", dest="verbose", help="Print little output"), 
    make_option("--voxs", type="character", default=NULL, help="A range of voxels to examine (this is mainly for testing purposes) and can be '1:10'."), 
    make_option("--ignoreprocerror", action="store_true", default=FALSE, help="Ignores the error generated if you specify the # of forks/threads to be greater than the actual number of estimated processes.")
)
    
# Make class/usage
parser <- OptionParser(usage = "%prog [options] output", option_list=option_list, add_help_option=TRUE)

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
  if (is.null(opts$jobs))
      set_parallel_procs(opts$forks, opts$threads, opts$verbose, opts$ignoreprocerror)
  # use foreach parallelization and shared memory?
  parallel_forks <- ifelse(opts$forks == 1, FALSE, TRUE)
  
  ###
  # Check Inputs
  ###
  vcat(opts$verbose, "Checking options")
  
  # check variables
  if (is.null(opts$indir))
      stop("must specify -i/--indir option")
  if (is.null(opts$model))
      stop("must specify -m/--model option")
  if (is.null(opts$formula))
      stop("must specify -f/--formula option")
  if (opts$overwrite)
      stop("Right now the overwrite function isn't implemented")
  if (opts$permutations < 2)
      stop("# of permutations must be greater than 1")
  
  # check paths exist
  ## indir
  opts$indir <- abspath(opts$indir)
  invisible(check_subdist(opts$indir))
  ## subdist
  if (is.null(opts$subdist))
      opts$subdist <- file.path(opts$indir, "subdist_gower.desc")
  if (!file.exists(opts$subdist))
      vstop("--subdist '%s' does not exist", opts$subdist)
  opts$subdist <- abspath(opts$subdist)
  ## model
  if (!file.exists(opts$model))
      stop("-m/--model ", opts$model, " does not exist")
  ## output
  opts$outdir <- file.path(opts$indir, args[1])
  if (file.exists(opts$outdir))
      stop("Output mdmr directory '", opts$outdir,  "' already exists")
  if (!file.exists(dirname(opts$outdir)))
      stop("Basepath to mdmr directory '", dirname(opts$outdir), "' doesn't exist")
  
  
  ###
  # Setup inputs
  ###
  vcat(opts$verbose, "Setting up inputs")
  
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
  if (is.null(opts$voxs)) {
      nvoxs <- ncol(tmp)
      voxs <- 1:nvoxs
  } else {
      voxs <- eval(parse(text=opts$voxs))
      nvoxs <- length(voxs)
  }
  rm(tmp); gc(FALSE, TRUE)
  
  # formula
  vcat(opts$verbose, "...formula")
  formula <- as.formula(sprintf("~ %s", opts$formula))
  vars <- as.character(as.list(attr(terms(formula), "variables"))[-1])
  
  # factors2perm
  vcat(opts$verbose, "...factors to permute")
  if (!is.null(opts$factors2perm)) {
      vcat(opts$verbose, "...factors2perm")
      
      opts$factors2perm <- sub(", ", ",", opts$factors2perm)
      opts$factors2perm <- strsplit(opts$factors2perm, ",")[[1]]
      
      for (x in opts$factors2perm) {
          if (!(x %in% vars)) {
              vstop("Factor to permute '%s' not found in formula (%s)", 
                    x, paste(vars, collapse=", "))
          }
      }
      
      nfactors <- length(opts$factors2perm)
  } else {
      nfactors <- length(attr(terms(formula), "term.labels"))
  }
  
  # model
  vcat(opts$verbose, "...model")
  model <- read.csv(opts$model)
  ## checks
  if (nrow(model) != nsubs)
      stop("# of rows in model file don't match # of subjects in distance matrix")
  for (v in vars) {
      if (is.null(model[[v]])) {
          vstop("Factor '%s' doesn't match any column in model file (%s)", 
                v, paste(colnames(model), collapse=", "))
      }
  }
  
  # strata
  if (!is.null(opts$strata)) {
      vcat(opts$verbose, "...strata")
      if (is.null(model[[opts$strata]]))
          stop("Strata given but doesn't match any column in the model file")
      else
          opts$strata <- model[[opts$strata]]
  }
    
  # output
  vcat(opts$verbose, "...creating output directory '%s'", opts$outdir)
  dir.create(opts$outdir)
  
  # subject distances
  vcat(opts$verbose, "Reading in subject distances")
  xdist <- attach.big.matrix(opts$subdist)
  xdist.path <- dirname(opts$subdist)
  
  # check
  vcat(opts$verbose, "...checking input")
  tmp <- matrix(xdist[,1], nsubs, nsubs)
  check_gmat(tmp)
  rm(tmp); invisible(gc(FALSE, TRUE))
  
  
  ###
  # Memory Demands: superblocksize & blocksize
  ###
  
  nperms <- opts$permutations + 1
  opts <- get_mdmr_memlimit(opts, nsubs, nvoxs, nperms, nfactors)
  
  
  ###
  # Compute MDMR
  ###
  start.time <- Sys.time()
  
  if (opts$saveperms) {
      fperms.path <- opts$outdir
  } else {
      fperms.path <- NULL
  }
  
  #save.image(file="z_env.rda")
  #q()
  
  if (is.null(opts$jobs)) {
      sge.info <- NULL
  } else {
      sge.info <- list(njobs=opts$jobs, nforks=opts$forks, nthreads=opts$threads, 
                       ignore.proc.error=opts$ignoreprocerror)
  }
  
  res.mdmr <- mdmr(xdist, formula, model, strata=opts$strata, 
                   nperms=opts$permutations, factors2perm=opts$factors2perm, 
                   superblocksize=opts$superblocksize, voxs=voxs, 
                   blocksize=opts$blocksize, verbose=verbosity, 
                   parallel=parallel_forks, shared=parallel_forks, 
                   G.path=xdist.path, fperms.path=fperms.path, save.fperms=TRUE, 
                   sge.info=sge.info)
  
  # Eventually remove calling of different functions for SGE vs not
  #if (is.null(opts$jobs)) {
  #    res.mdmr <- mdmr(xdist, formula, model, nperms=opts$permutations, 
  #                     superblocksize=opts$superblocksize, blocksize=opts$blocksize, 
  #                     strata=opts$strata, factors2perm=opts$factors2perm, 
  #                     verbose=verbosity, parallel=parallel_forks, 
  #                     G.path=xdist.path, fperms.path=fperms.path, 
  #                     voxs=voxs)
  #} else {
  #    res.mdmr <- mdmr.sge(opts$subdist, formula, model, nperms=opts$permutations, 
  #                     superblocksize=opts$superblocksize, blocksize=opts$blocksize, 
  #                     strata=opts$strata, factors2perm=opts$factors2perm, 
  #                     verbose=verbosity, parallel=parallel_forks, 
  #                     fperms.path=fperms.path, 
  #                     forks=opts$forks, threads=opts$threads, njobs=opts$jobs, 
  #                     voxs=voxs, ignore.proc.error=opts$ignoreprocerror)
  #}
  rm(xdist)
  invisible(gc(FALSE, TRUE))
  
  end.time <- Sys.time()
  vcat(opts$verbose, "MDMR is done! It took: %.2f minutes\n", 
       as.numeric(end.time-start.time, units="mins"))
       
   #options(error=recover)
   save(res.mdmr, file="z_res_mdmr.rda")
   #res.mdmr$pvals
   #save(res.mdmr$pvals, file="ztest2.rda")
  
  
  ###
  # Save MDMR Results
  ###
  save_mdmr(res.mdmr, voxs, opts$indir, opts$outdir, opts$verbose)
  rm(res.mdmr)
  invisible(gc(FALSE, TRUE))

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
