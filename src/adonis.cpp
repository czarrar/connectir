#include "connectir/connectir.h"

SEXP mdmr_fstats_to_pvals(SEXP SFmat)
{
    try {
        arma::mat Fmat(1,1);
        const double* old_fptr = sbm_to_arma_xd(SFmat, Fmat);
        
        double nvoxs = static_cast<double>(Fmat.n_cols);
        double nperms = static_cast<double>(Fmat.n_rows);
        
        Rcpp::NumericVector pvals(nvoxs);
        
        // original F-stats
        arma::rowvec realFs = Fmat.row(0);
        
        double i;
        for (i = 0; i < nvoxs; ++i)
        {
            pvals[i] = arma::as_scalar(arma::sum(Fmat.col(i) >= realFs(i))/nperms);
        }
        
        free_arma(Fmat, old_fptr);
        
        return Rcpp::wrap( pvals );;
    } catch(std::exception &ex) {
        forward_exception_to_r(ex);
    } catch(...) {
        ::Rf_error("c++ exception (unknown reason)");
    }
    
    return R_NilValue;
}

SEXP mdmr_worker(SEXP SGmat, SEXP SFperms, 
                 SEXP SH2mats, SEXP SIHmats, 
                 SEXP SdfRes, SEXP SdfExp)
{
    try {
        // nterms
        Rcpp::List tmp(SH2mats);
        index_type nterms = static_cast<index_type>(tmp.size());
        
        // Gs and nvoxs
        arma::mat Gmat(1,1);
        const double* old_gptr = sbm_to_arma_xd(SGmat, Gmat);
        index_type nvoxs = static_cast<index_type>(Gmat.n_cols);

        // Fperms
        SEXP SFmat;
        arma::mat Fmat(1,1); const double* old_fptr = Fmat.memptr();
        
        // H2s
        SEXP SH2mat;
        arma::mat H2mat(1,1); const double* old_h2ptr = H2mat.memptr();
        
        // IHs
        SEXP SIHmat;
        arma::mat IHmat(1,1); const double* old_ihptr = IHmat.memptr();
        
        // dfs
        double dfRes = DOUBLE_DATA(SdfRes)[0];
        Rcpp::NumericVector RdfExp(SdfExp);
        arma::vec dfExp(RdfExp.begin(), RdfExp.size(), false);
        
        arma::mat ExplainedVariance;
        arma::mat ErrorVariance;
        index_type i, j;
        for (i=0; i < nterms; ++i)
        {
            
            /***
            * Explained Variance
            ***/
            
            // H2
            PROTECT(SH2mat = VECTOR_ELT(SH2mats, i));
            sbm_to_arma_xd(SH2mat, H2mat);
            UNPROTECT(1);
            
            // Explained Variance
            ExplainedVariance = arma::trans(H2mat) * Gmat;
            
            
            /***
            * Error Variance
            ***/
            
            // IH
            PROTECT(SIHmat = VECTOR_ELT(SIHmats, i));
            sbm_to_arma_xd(SIHmat, IHmat);
            UNPROTECT(1);
            
            // Error Variance
            arma::mat ErrorVariance = arma::trans(IHmat) * Gmat;
            
            /***
            * Peusod-F Statistic
            ***/
            
            // Fstats (output)
            PROTECT(SFmat = VECTOR_ELT(SFperms, i));
            sbm_to_arma_xd(SFmat, Fmat);
            UNPROTECT(1);
            
            // ExplainedVariance/ErrorVariance
            Fmat = (ExplainedVariance/ErrorVariance) * (dfRes/dfExp(i));
            
        }
        
        // Clear matrices
        free_arma(Gmat, old_gptr);
        free_arma(H2mat, old_h2ptr);
        free_arma(IHmat, old_ihptr);
        free_arma(Fmat, old_fptr);
        
        return R_NilValue;
    } catch(std::exception &ex) {
        forward_exception_to_r(ex);
    } catch(...) {
        ::Rf_error("c++ exception (unknown reason)");
    }
    
    return R_NilValue;
}

template<typename CType, typename BMAccessorType>
SEXP ComputePvals(BigMatrix *inMat, BigMatrix *outMat, double colnum) {
    BMAccessorType om( *outMat );
    CType *outCol;
    outCol = om[static_cast<index_type>(colnum)-1];
    
    BMAccessorType im( *inMat );    
    index_type ncols = inMat->ncol();
    index_type nrows = inMat->nrow();
    index_type i=0;
    index_type j=0;
    CType *inCol;
    
    for (i=0; i < ncols; ++i) {
        inCol = im[i];
        outCol[i] = 1;
        for (j=1; j < nrows; ++j) {
            if (inCol[j] > inCol[0])
                outCol[i] += 1;
        }
        outCol[i] = outCol[i]/nrows;
    }
    
    return R_NilValue;
}


extern "C" {

SEXP ComputePvalsMain(SEXP Rinmat, SEXP Routmat, SEXP Routcol) {
    BigMatrix *inMat = reinterpret_cast<BigMatrix*>(R_ExternalPtrAddr(Rinmat));
    BigMatrix *outMat = reinterpret_cast<BigMatrix*>(R_ExternalPtrAddr(Routmat));
    double outCol = NUMERIC_DATA(Routcol)[0];
    
    if (inMat->separated_columns() != outMat->separated_columns())
        Rf_error("all big matrices are not the same column separated type");
    if (inMat->matrix_type() != outMat->matrix_type())
        Rf_error("all big matrices are not the same matrix type");
    if (inMat->ncol() != outMat->nrow())
        Rf_error("inMat # of cols must be the same as outMat # of rows");
    
    CALL_BIGFUNCTION_ARGS_THREE(ComputePvals, inMat, outMat, outCol)
    return(ret);
}

} // end extern "C"
