library(RcppArmadillo)
library(inline)

plugin_bigmemory <- function() {
  l <- getPlugin("RcppArmadillo")
  
  l$includes <- paste(l$includes, '
#include "bigmemory/BigMatrix.h"
#include "bigmemory/MatrixAccessor.hpp"
#include "bigmemory/bigmemoryDefines.h"
#include "bigmemory/isna.hpp"
')
  
  l$LinkingTo <- c("bigmemory", l$LinkingTo)
  
  l$Depends <- c("bigmemory", l$Depends)  
  
  return(l)
}
registerPlugin("bigmemory", plugin_bigmemory)

plugin_connectir <- function() {
  l <- getPlugin("bigmemory")
  
  l$includes <- paste(l$includes, '
#include "connectir/connectir.h"
')
  
  l$LinkingTo <- c("connectir", l$LinkingTo)
  
  l$Depends <- c("connectir", l$Depends)  
  
  return(l)
}
registerPlugin("connectir", connectir)


