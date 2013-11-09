#include "connectir/connectir.h"

// Get rank of big matrix
SEXP big_qlm_rank(SEXP Xr) {    
    try {
        BM_TO_ARMA_ONCE(Xr, X)
        using namespace arma;
        return Rcpp::wrap(rank(X));
    } catch(std::exception &ex) {
        forward_exception_to_r(ex);
    } catch(...) {
        ::Rf_error("c++ exception (unknown reason)");
    }    
    return R_NilValue;
}

// Get diag( inv(t(X) %*% X) )
SEXP big_qlm_dd(SEXP Xr) {    
    try {
        BM_TO_ARMA_ONCE(Xr, X)
        using namespace arma;
        return Rcpp::wrap( diagvec( inv(trans(X) * X) ) );
    } catch(std::exception &ex) {
        forward_exception_to_r(ex);
    } catch(...) {
        ::Rf_error("c++ exception (unknown reason)");
    }    
    return R_NilValue;
}

// Solve for y ~ X (also pass outputs here)
SEXP big_qlm_fit(SEXP yr, SEXP Xr, SEXP coefr, SEXP residr, SEXP mser, SEXP nr, SEXP kr, SEXP mr) { 
    try {
        BM_TO_ARMA_INIT()
        BM_TO_ARMA_MULTIPLE(yr, y)
        BM_TO_ARMA_MULTIPLE(Xr, X)
        BM_TO_ARMA_MULTIPLE(coefr, coef)
        BM_TO_ARMA_MULTIPLE(residr, resid)
        BM_TO_ARMA_MULTIPLE(mser, mse)
        
        //NumericVector nc(nr);
        //NumericVector kc(kr);
        //NumericVector mc(mr);
        double n = DOUBLE_DATA(nr)[0];
        double k = DOUBLE_DATA(kr)[0];
        double m = DOUBLE_DATA(mr)[0];
        
        // fit model y ~ X
        coef = arma::solve(X, y);
        
        // residuals
        resid = y - X*coef;
        
        // mean-square of errors
        mse = arma::sum( arma::pow(resid, 2), 0 )/(n-k);
    } catch(std::exception &ex) {
        forward_exception_to_r(ex);
    } catch(...) {
        ::Rf_error("c++ exception (unknown reason)");
    }    
    return R_NilValue;
}

SEXP big_qlm_rsquared(SEXP yr, SEXP coefr, SEXP residr, SEXP mser, SEXP nr, SEXP kr, 
                      SEXP rsquaredr, SEXP to_adjustr) 
{
    try {
        BM_TO_ARMA_INIT()
        BM_TO_ARMA_MULTIPLE(yr, y)
        BM_TO_ARMA_MULTIPLE(coefr, coef)
        BM_TO_ARMA_MULTIPLE(residr, resid)
        BM_TO_ARMA_MULTIPLE(mser, mse)
        BM_TO_ARMA_MULTIPLE(rsquaredr, rsquared)
        
        double n = DOUBLE_DATA(nr)[0];
        double k = DOUBLE_DATA(kr)[0];
        int to_adjust = INTEGER_DATA(to_adjustr)[0];
        
        // fitted data & mean sum squares
        arma::mat fitted = y - resid;
        
        arma::rowvec mfit = arma::mean(fitted, 0);
        
        for(size_t i = 0; i < fitted.n_cols; ++i)
            fitted.col(i) = fitted.col(i) - mfit(i);
        
        arma::rowvec mss = arma::sum(arma::pow( fitted, 2), 0);
        
        // residual sum squares
        arma::rowvec rss = arma::trans( mse * (n-k) );
        
        // r2
        rsquared = mss/(mss+rss);
        
        // adjusted r2
        if (to_adjust == 1)
            rsquared = 1 - (1 - rsquared) * ((n-1)/(n-k));
    } catch(std::exception &ex) {
        forward_exception_to_r(ex);
    } catch(...) {
        ::Rf_error("c++ exception (unknown reason)");
    }    
    return R_NilValue;
}

// Solve for y ~ X (also pass outputs here)
SEXP big_qlm_residuals(SEXP yr, SEXP Xr, SEXP residr, SEXP add_meanr) { 
    try {
        BM_TO_ARMA_INIT()
        BM_TO_ARMA_MULTIPLE(yr, y)
        BM_TO_ARMA_MULTIPLE(Xr, X)
        BM_TO_ARMA_MULTIPLE(residr, resid)
			
		int add_mean = INTEGER_DATA(add_meanr)[0];
        
        // fit model y ~ X
        arma::mat coef = arma::solve(X, y);
        
        // residuals
        resid = y - X*coef;
		
		// add back mean
		if (add_mean) {
			arma::rowvec m = mean(X);
	        for(size_t i = 0; i < X.n_cols; ++i)
				resid.col(i) = resid.col(i) + m(i);
		}
        
        return R_NilValue;
    } catch(std::exception &ex) {
        forward_exception_to_r(ex);
    } catch(...) {
        ::Rf_error("c++ exception (unknown reason)");
    }    
    return R_NilValue;
}

// Solve for y ~ X (also pass outputs here)
SEXP big_qlm_contrasts(SEXP fit_coefr, SEXP fit_mser, SEXP conr, SEXP ddr, SEXP coefr, SEXP ser, SEXP tvalr, SEXP mr) { 
    try {
        BM_TO_ARMA_INIT()
        BM_TO_ARMA_MULTIPLE(fit_coefr, fit_coef)
        BM_TO_ARMA_MULTIPLE(fit_mser, fit_mse)
        BM_TO_ARMA_MULTIPLE(coefr, coef)
        BM_TO_ARMA_MULTIPLE(ser, se)
        BM_TO_ARMA_MULTIPLE(tvalr, tval)
        
        Rcpp::NumericMatrix conc(conr);
        arma::mat con(conc.begin(), conc.nrow(), conc.ncol(), false);
        Rcpp::NumericVector ddc(ddr);
        arma::colvec dd(ddc.begin(), ddc.size(), false);
        
        double m = DOUBLE_DATA(mr)[0];
        
        coef = con * fit_coef;
        
        arma::mat c_dd = con * dd;
        for(size_t i = 0; i < m; ++i)
        {
          se.col(i) = arma::sqrt( fit_mse[i] * c_dd);
        }
        
        tval = coef/se;
    } catch(std::exception &ex) {
        forward_exception_to_r(ex);
    } catch(...) {
        ::Rf_error("c++ exception (unknown reason)");
    }    
    return R_NilValue;
}


