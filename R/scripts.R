# Class/Functions that are meant to be called by the connectir scripts

CONNECTIR <- setRefClass("connectir",
    fields = list(
        opts = "list"
    )
)


###
# Misc
###

printf <- function(msg, ..., newline=TRUE) {
    if (opts$verbose) {
        cat(sprintf(msg, ...))
        if (newline) cat("\n")
    }
}


###
# Parallel Processing
###

parallelize <- function() {
    printf("Setting %i cores to be used", opts$cores)
    if (opts$cores > 1) {
        printf(opts, "...setting parallel processing with doMC")
        suppressPackageStartupMessages(library("doMC"))
        registerDoMC()
        if (opts$cores > getDoParWorkers())
        	stop("Number of -c/--cores specified '", opts$cores, "' is greater than the actual number of cores '", getDoParWorkers(), "'")
    }
    options(cores=opts$cores)
    
    printf("04. Setting %i MKL threads to be used", opts$threads)
    if (length(grep("bigalgebra", search())) > 0)
        detach(pos=grep("bigalgebra", search()), unload=TRUE)
    if (existsFunction("setMKLthreads"))
    	setMKLthreads(opts$threads)
    else
        Sys.setenv(MKL_NUM_THREADS=opts$threads)
    suppressPackageStartupMessages(library("bigalgebra"))
}


