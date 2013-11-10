#include "connectir/connectir.h"

// Return row sum of big matrix (must be of type double)
SEXP big_rowsum(SEXP SX, SEXP SX_firstCol, SEXP SX_lastCol) {
    try {
        arma::mat X;
        sub_sbm_to_arma_xd(SX, X, SX_firstCol, SX_lastCol);
        return Rcpp::wrap(sum(X, 1) );
    } catch(std::exception &ex) {
        forward_exception_to_r(ex);
    } catch(...) {
        ::Rf_error("c++ exception (unknown reason)");
    }
    return R_NilValue;
}

// Return row sum of big matrix (must be of type double)
SEXP big_rowmean(SEXP SX, SEXP SX_firstCol, SEXP SX_lastCol) {
    try {
        arma::mat X;
        sub_sbm_to_arma_xd(SX, X, SX_firstCol, SX_lastCol);
        return Rcpp::wrap(mean(X, 1) );
    } catch(std::exception &ex) {
        forward_exception_to_r(ex);
    } catch(...) {
        ::Rf_error("c++ exception (unknown reason)");
    }
    return R_NilValue;
}

// Return row sum of big matrix (must be of type double)
SEXP big_colsum(SEXP SX, SEXP SX_firstCol, SEXP SX_lastCol) {
    try {
        arma::mat X;
        sub_sbm_to_arma_xd(SX, X, SX_firstCol, SX_lastCol);
        return Rcpp::wrap(sum(X, 0) );
    } catch(std::exception &ex) {
        forward_exception_to_r(ex);
    } catch(...) {
        ::Rf_error("c++ exception (unknown reason)");
    }
    return R_NilValue;
}

// Return row sum of big matrix (must be of type double)
SEXP big_colmean(SEXP SX, SEXP SX_firstCol, SEXP SX_lastCol) {
    try {
        arma::mat X;
        sub_sbm_to_arma_xd(SX, X, SX_firstCol, SX_lastCol);
        return Rcpp::wrap(mean(X, 0) );
    } catch(std::exception &ex) {
        forward_exception_to_r(ex);
    } catch(...) {
        ::Rf_error("c++ exception (unknown reason)");
    }
    return R_NilValue;
}
