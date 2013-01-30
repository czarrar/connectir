#include "connectir/connectir.h"

SEXP subdist_combine_and_scale_submaps(SEXP Slist_corMaps, SEXP Sseed, SEXP SvoxInds, 
                                       SEXP SseedCorMaps)
{
    try {
        // Setup Inputs
        Rcpp::List list_corMaps(Slist_corMaps);
        index_type seed = static_cast<index_type>(DOUBLE_DATA(Sseed)[0] - 1);
        Rcpp::NumericVector voxInds(SvoxInds);
        index_type nsubs = static_cast<index_type>(list_corMaps.size());
        index_type nvoxs = static_cast<index_type>(voxInds.size());
        arma::mat seedCorMaps(1,1); const double* old_mptr = seedCorMaps.memptr();
        sbm_to_arma_xd(SseedCorMaps, seedCorMaps); double* seedMap;
        
        // Copy over subject connectivity maps from list to matrix
        // and scale
        SEXP SsubCorMaps; 
        arma::mat subCorMaps(1,1); const double* old_sptr = subCorMaps.memptr();
        index_type voxi, sub, vox;
        LDOUBLE x, delta, mean, M2, stdev, scaled_x;
        for (sub = 0; sub < nsubs; ++sub)
        {
            PROTECT(SsubCorMaps = VECTOR_ELT(Slist_corMaps, sub));
            sbm_to_arma_xd(SsubCorMaps, subCorMaps);
            UNPROTECT(1);
            
            delta = mean = M2 = stdev = 0;
            for (voxi=0; voxi < nvoxs; ++voxi) {
                // todo: add checking for NaN...but shouldn't really have any!
                vox = static_cast<index_type>(voxInds(voxi)-1);
                x = static_cast<LDOUBLE>(subCorMaps(seed,vox));
                delta = x - mean;
                mean = mean + delta/static_cast<LDOUBLE>(voxi+1);
                M2 = M2 + delta*(x - mean);
            }
            stdev = sqrt(M2/(static_cast<LDOUBLE>(nvoxs-1)));
            
            seedMap = const_cast<double *>(seedCorMaps.colptr(sub));
            for (voxi=0; voxi < nvoxs; ++voxi) {
                vox = static_cast<index_type>(voxInds(voxi)-1);
                scaled_x = (static_cast<LDOUBLE>(subCorMaps(seed,vox))-mean)/stdev;
                seedMap[voxi] = static_cast<double>(scaled_x);
            }
        }
        
        seedMap = NULL;
        free_arma(seedCorMaps, old_mptr);
        free_arma(subCorMaps, old_sptr);
        
        return SseedCorMaps;
    } catch(std::exception &ex) {
        forward_exception_to_r(ex);
    } catch(...) {
        ::Rf_error("c++ exception (unknown reason)");
    }
    
    return R_NilValue;
}

SEXP subdist_combine_submaps(SEXP Slist_corMaps, SEXP Sseed, SEXP SvoxInds, 
                             SEXP SseedCorMaps)
{
    try {
        // Setup Inputs
        Rcpp::List list_corMaps(Slist_corMaps);
        index_type seed = static_cast<index_type>(DOUBLE_DATA(Sseed)[0] - 1);
        Rcpp::NumericVector voxInds(SvoxInds);
        index_type nsubs = static_cast<index_type>(list_corMaps.size());
        index_type nvoxs = static_cast<index_type>(voxInds.size());
        arma::mat seedCorMaps(1,1); const double* old_mptr = seedCorMaps.memptr();
        sbm_to_arma_xd(SseedCorMaps, seedCorMaps); double* seedMap;
        
        // Copy over subject connectivity maps from list to matrix
        // and scale
        SEXP SsubCorMaps; 
        arma::mat subCorMaps(1,1); const double* old_sptr = subCorMaps.memptr();
        index_type voxi, sub, vox;
        for (sub = 0; sub < nsubs; ++sub)
        {
            PROTECT(SsubCorMaps = VECTOR_ELT(Slist_corMaps, sub));
            sbm_to_arma_xd(SsubCorMaps, subCorMaps);
            UNPROTECT(1);
            
            seedMap = const_cast<double *>(seedCorMaps.colptr(sub));
            for (voxi=0; voxi < nvoxs; ++voxi) {
                vox = static_cast<index_type>(voxInds(voxi)-1);
                seedMap[voxi] = subCorMaps(seed,vox);
            }
        }
        
        seedMap = NULL;
        free_arma(seedCorMaps, old_mptr);
        free_arma(subCorMaps, old_sptr);
        
        return SseedCorMaps;
    } catch(std::exception &ex) {
        forward_exception_to_r(ex);
    } catch(...) {
        ::Rf_error("c++ exception (unknown reason)");
    }
    
    return R_NilValue;
}

SEXP subdist_combine_and_trans_submaps(SEXP Slist_corMaps, SEXP Sseed, SEXP SvoxInds, 
                                       SEXP SseedCorMaps)
{
    try {
        // Setup Inputs
        Rcpp::List list_corMaps(Slist_corMaps);
        index_type seed = static_cast<index_type>(DOUBLE_DATA(Sseed)[0] - 1);
        Rcpp::NumericVector voxInds(SvoxInds);
        index_type nsubs = static_cast<index_type>(list_corMaps.size());
        index_type nvoxs = static_cast<index_type>(voxInds.size());
        arma::mat seedCorMaps(1,1); const double* old_mptr = seedCorMaps.memptr();
        sbm_to_arma_xd(SseedCorMaps, seedCorMaps);
        
        // Copy over subject connectivity maps from list to matrix
        // and scale
        SEXP SsubCorMaps; 
        arma::mat subCorMaps(1,1); const double* old_sptr = subCorMaps.memptr();
        index_type voxi, sub, vox;
        for (sub = 0; sub < nsubs; ++sub)
        {
            PROTECT(SsubCorMaps = VECTOR_ELT(Slist_corMaps, sub));
            sbm_to_arma_xd(SsubCorMaps, subCorMaps);
            UNPROTECT(1);
            
            for (voxi=0; voxi < nvoxs; ++voxi) {
                vox = static_cast<index_type>(voxInds(voxi)-1);
                seedCorMaps(sub,voxi) = subCorMaps(seed,vox);
            }
        }
        
        free_arma(seedCorMaps, old_mptr);
        free_arma(subCorMaps, old_sptr);
        
        return SseedCorMaps;
    } catch(std::exception &ex) {
        forward_exception_to_r(ex);
    } catch(...) {
        ::Rf_error("c++ exception (unknown reason)");
    }
    
    return R_NilValue;
}


// 1 - pearson correlation
SEXP subdist_pearson_distance(SEXP SseedCorMaps, SEXP Sdmats, SEXP Sdcol, 
                              SEXP Sistrans)
{
    try {
        arma::mat seedCorMaps(1,1);
        const double* old_sptr = sbm_to_arma_xd(SseedCorMaps, seedCorMaps);
                
        arma::mat dmat(1,1);
        const double* old_dptr = sub_sbm_to_arma_xd(Sdmats, dmat, Sdcol, Sdcol);
        index_type nsubs = static_cast<index_type>(
                                sqrt(static_cast<double>(dmat.n_rows)));
        dmat.reshape(nsubs, nsubs);
        
        Rcpp::LogicalVector istrans(Sistrans);
        if (istrans[0] == TRUE) {
            index_type nvoxs = static_cast<index_type>(seedCorMaps.n_cols);
            index_type df = nvoxs - 1;
            dmat = 1 - (seedCorMaps * arma::trans(seedCorMaps))/df;
        } else {
            index_type nvoxs = static_cast<index_type>(seedCorMaps.n_rows);
            index_type df = nvoxs - 1;
            dmat = 1 - (arma::trans(seedCorMaps) * seedCorMaps)/df;
        }
        
        free_arma(seedCorMaps, old_sptr);
        free_arma(dmat, old_dptr);
        
        return Sdmats;
    } catch(std::exception &ex) {
        forward_exception_to_r(ex);
    } catch(...) {
        ::Rf_error("c++ exception (unknown reason)");
    }
    
    return R_NilValue;
}

// Y <- -0.5 * dmat^2 %*% (I - ones %*% t(ones)/n)
SEXP big_gower(SEXP SX, SEXP SY, 
               SEXP SX_firstCol, SEXP SX_lastCol, 
               SEXP SY_firstCol, SEXP SY_lastCol)
{
    try {
        arma::mat X(1,1);
        const double* old_Xptr = sub_sbm_to_arma_xd(SX, X, SX_firstCol, SX_lastCol);
        arma::mat Y(1,1);
        const double* old_Yptr = sub_sbm_to_arma_xd(SY, Y, SY_firstCol, SY_lastCol);
                
        if (X.n_rows != Y.n_rows || X.n_cols != Y.n_cols)
            Rf_error("dimension mismatch between input and output matrices");
        index_type n = static_cast<index_type>(sqrt(
                        static_cast<double>(Y.n_rows)));
        
        using namespace arma;
        
        // adj = I - ones %*% t(ones)/n
        mat I = eye<mat>(n, n);
        colvec z_ones = ones<colvec>(n);
        mat adj = I - z_ones * trans(z_ones) / n;
        
        // A = -0.5 * dmat^2
        Y = arma::pow(X, 2)/(-2);
        
        // G = A %*% adj
        mat Y_vox(n,n);
        const double *old_ptr = Y_vox.memptr();
        for (index_type i = 0; i < Y.n_cols; ++i)
        {
            arma::access::rw(Y_vox.mem) = Y.colptr(i);
            Y_vox = Y_vox * adj;
        }
        
        free_arma(X, old_Xptr);
        free_arma(Y, old_Yptr);
        arma::access::rw(Y_vox.mem) = old_ptr;
        
        return Rcpp::wrap( SY ); // or return R_NilValue?
    } catch(std::exception &ex) {
        forward_exception_to_r(ex);
    } catch(...) {
        ::Rf_error("c++ exception (unknown reason)");
    }
    return R_NilValue;
}
