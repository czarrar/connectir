# qlm => quick multiple linear regression

# Script that will go through each connectivity map and compute regression analysis
#library(inline)
#
#inline_qlm_rank <- cxxfunction( signature(Xs = "numeric"), 
#'
#  using namespace Rcpp;
#  
#  try{
#    NumericMatrix Xr(Xs); // design matrix
#    arma::mat X(Xr.begin(), Xr.nrow(), Xr.ncol(), false);
#    return wrap( arma::rank( X ) );
#  } catch( std::exception &ex ) {
#    forward_exception_to_r( ex );
#  } catch(...) { 
#    ::Rf_error( "c++ exception (unknown reason)" ); 
#  }
#  return R_NilValue; // -Wall
#', plugin = "RcppArmadillo")
#
#inline_qlm_dd <- cxxfunction( signature(Xs = "numeric"), 
#'
#  using namespace Rcpp;
#  
#  try{
#    NumericMatrix Xr(Xs); // design matrix
#    arma::mat X(Xr.begin(), Xr.nrow(), Xr.ncol(), false);
#    return wrap( arma::diagvec( arma::inv(arma::trans(X) * X) ) );
#  } catch( std::exception &ex ) {
#    forward_exception_to_r( ex );
#  } catch(...) { 
#    ::Rf_error( "c++ exception (unknown reason)" ); 
#  }
#  return R_NilValue; // -Wall
#', plugin = "RcppArmadillo")
#
#inline_qlm_fit <- cxxfunction( signature(ys = "numeric", Xs = "numeric"), 
#'
#  using namespace Rcpp;
#
#  try {
#    NumericMatrix yr(ys);   // data
#    NumericMatrix Xr(Xs);   // design matrix
#    
#    int n = Xr.nrow(), k = Xr.ncol();
#    int m = yr.ncol();
#    if (n != yr.nrow())
#      ::Rf_error("size mismatch. X must have the same number of row as y (where y ~ X).");
#    
#    arma::mat X(Xr.begin(), n, k, false);
#    arma::mat y(yr.begin(), n, m, false);
#    
#    // fit model y ~ X
#    arma::mat coefs = arma::solve(X, y);
#    
#    // residuals
#    arma::mat resids = y - X*coefs;
#    
#    // mean-square of errors
#    arma::rowvec mse = arma::sum( arma::pow(resids, 2), 0 )/(n-k);
#    
#    // TODO: have option where can calculate standard error + tvals
#    
#    return List::create(
#        _["coefficients"] = coefs,
#        _["residuals"]    = resids,
#        _["mse"]          = mse
#    );
#  } catch( std::exception &ex ) {
#  forward_exception_to_r( ex );
#  } catch(...) { 
#  ::Rf_error( "c++ exception (unknown reason)" ); 
#  }
#  return R_NilValue; // -Wall
#', plugin = "RcppArmadillo")
#
#inline_qlm_contrasts <- cxxfunction( signature( coefs_s="numeric", mse_s="numeric", dd_s="numeric", cons_s="numeric" ), 
#'
#  using namespace Rcpp;
#
#  try{
#    NumericMatrix coefs_r(coefs_s);
#    NumericVector mse_r(mse_s);
#    NumericVector dd_r(dd_s);
#    NumericMatrix cons_r(cons_s);
#  
#    int ne = coefs_r.nrow();  // # of regressors
#    int ni = coefs_r.ncol();  // # of ys (i.e., # of independent regression analyses)
#    int nc = cons_r.nrow();    // # of contrasts
#    if (ne != cons_r.ncol())
#      ::Rf_error("size mismatch between Coefficients and Contrasts.");
#    if (ni != mse_r.size())
#      ::Rf_error("size mismatch between Coefficients and MSE.");
#    
#    arma::mat coefs(coefs_r.begin(), ne, ni, false);
#    arma::rowvec mse(mse_r.begin(), ni, false);
#    arma::colvec dd(dd_r.begin(), dd_r.size(), false);
#    arma::mat cons(cons_r.begin(), nc, ne, false);
#  
#    arma::mat contrast_coefs = cons * coefs;
#    arma::mat contrast_dd = cons * dd;
#  
#    arma::mat contrast_se(nc, ni);
#    for(size_t i = 0; i < ni; ++i)
#    {
#      contrast_se.col(i) = arma::sqrt( mse[i] * contrast_dd);
#    }
#  
#    arma::mat tvals = contrast_coefs/contrast_se;
#  
#    return List::create(
#        _["coefficients"]     = contrast_coefs, 
#        _["standard_errors"]  = contrast_se, 
#        _["tvals"]            = tvals
#    );
#  } catch( std::exception &ex ) {
#    forward_exception_to_r( ex );
#  } catch(...) { 
#    ::Rf_error( "c++ exception (unknown reason)" ); 
#  }
#  return R_NilValue;
#', plugin="RcppArmadillo")
#
#
#