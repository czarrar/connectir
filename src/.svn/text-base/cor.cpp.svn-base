#include "connectir/connectir.h"

// In R:
// library(connectir)
// x <- as.big.matrix(matrix(1:20, 5, 4), type="double")
// y <- .Call("test_sub_matrix", x, as.double(2), as.double(4))
// all.equal(x[,2:4], y)
SEXP test_sub_matrix(SEXP As, SEXP As_firstCol, SEXP As_lastCol) {
    try {
        BM_COL_INIT(As_firstCol, A_firstCol)
        BM_COL_INIT(As_lastCol, A_lastCol)
    
        BM_TO_ARMA_INIT()
        
        SUB_BM_TO_ARMA_MULTIPLE(As, A, A_firstCol, A_lastCol)
        
        return Rcpp::wrap( A );
    } catch(std::exception &ex) {
        forward_exception_to_r(ex);
    } catch(...) {
        ::Rf_error("c++ exception (unknown reason)");
    }
    return R_NilValue;
}

// Cs <- (t(As) %*% Bs)/(nrow(Bs)-1)
SEXP big_cor(SEXP As, SEXP Bs, SEXP Cs, 
             SEXP As_firstCol, SEXP As_lastCol, 
             SEXP Bs_firstCol, SEXP Bs_lastCol, 
             SEXP Cs_firstCol, SEXP Cs_lastCol) 
{    
    try {
        BM_COL_INIT(As_firstCol, A_firstCol)
        BM_COL_INIT(As_lastCol, A_lastCol)
        BM_COL_INIT(Bs_firstCol, B_firstCol)
        BM_COL_INIT(Bs_lastCol, B_lastCol)
        BM_COL_INIT(Cs_firstCol, C_firstCol)
        BM_COL_INIT(Cs_lastCol, C_lastCol)
        
        BM_TO_ARMA_INIT()
        SUB_BM_TO_ARMA_MULTIPLE(As, A, A_firstCol, A_lastCol)
        SUB_BM_TO_ARMA_MULTIPLE(Bs, B, B_firstCol, B_lastCol)
        double df = pMat->nrow() - 1;
        SUB_BM_TO_ARMA_MULTIPLE(Cs, C, C_firstCol, C_lastCol)
        
        // TODO: make this an explicit user option!
        if (ncol == 1 && B.n_cols != 1 && C.n_rows == (A.n_cols*B.n_cols))
            C.reshape(A.n_cols, B.n_cols);
        
        if (A.n_cols != C.n_rows)
            Rf_error("incorrect number of rows in output matrix");
        if (B.n_cols != C.n_cols)
            Rf_error("incorrect number of columns in output matrix");
        
        C = (arma::trans(A) * B)/df;
        
        return R_NilValue;
    } catch(std::exception &ex) {
        forward_exception_to_r(ex);
    } catch(...) {
        ::Rf_error("c++ exception (unknown reason)");
    }
    return R_NilValue;
}

// Cs <- (As %*% t(Bs))/(nrow(As)-1)
SEXP big_tcor(SEXP As, SEXP Bs, SEXP Cs, 
             SEXP As_firstCol, SEXP As_lastCol, 
             SEXP Bs_firstCol, SEXP Bs_lastCol, 
             SEXP Cs_firstCol, SEXP Cs_lastCol) 
{    
    try {
        BM_COL_INIT(As_firstCol, A_firstCol)
        BM_COL_INIT(As_lastCol, A_lastCol)
        BM_COL_INIT(Bs_firstCol, B_firstCol)
        BM_COL_INIT(Bs_lastCol, B_lastCol)
        BM_COL_INIT(Cs_firstCol, C_firstCol)
        BM_COL_INIT(Cs_lastCol, C_lastCol)
        
        BM_TO_ARMA_INIT()
        SUB_BM_TO_ARMA_MULTIPLE(As, A, A_firstCol, A_lastCol)
        SUB_BM_TO_ARMA_MULTIPLE(Bs, B, B_firstCol, B_lastCol)
        double df = ncol - 1;
        SUB_BM_TO_ARMA_MULTIPLE(Cs, C, C_firstCol, C_lastCol)
        
        // TODO: make this an explicit user option!
        if (ncol == 1 && B.n_rows != 1 && C.n_rows == (A.n_rows*B.n_rows))
            C.reshape(A.n_rows, B.n_rows);
        
        if (A.n_rows != C.n_rows)
            Rf_error("incorrect number of rows in output matrix");
        if (B.n_rows != C.n_cols)
            Rf_error("incorrect number of columns in output matrix");
        
        C = (A * arma::trans(B))/df;
        
        return R_NilValue;
    } catch(std::exception &ex) {
        forward_exception_to_r(ex);
    } catch(...) {
        ::Rf_error("c++ exception (unknown reason)");
    }
    
    return R_NilValue;
}

double kendall_worker_d(SEXP SRatings)
{
    try {
        arma::mat Ratings(1,1);
        const double* old_rptr = sbm_to_arma_xd(SRatings, Ratings);
    
        arma::uvec Ovec(Ratings.n_rows);
        arma::mat Rmat(Ratings.n_rows, Ratings.n_cols);
        double* Rvec;
    
        // Ranking
        arma::u32 i, j;
        for (i = 0; i < Ratings.n_cols; ++i) {
            Ovec = arma::sort_index(Ratings.unsafe_col(i));
            Rvec = Rmat.colptr(i);
            for (j = 0; j < Ratings.n_rows; ++j) {
                Rvec[Ovec(j)] = static_cast<double>(j+1);
            }
        }
    
        // rowSums
        arma::colvec rs = arma::sum(Rmat, 1);
    
        // var
        double v = var(rs);
    
        // clean up
        Rvec = NULL;
        free_arma(Ratings, old_rptr);
    
        return v;
    } catch( std::exception &ex ) {
      forward_exception_to_r( ex );
    } catch(...) { 
      ::Rf_error( "c++ exception (unknown reason)" ); 
    }
    return -1;
}

SEXP kendall_worker(SEXP SRatings)
{
    try {
        return Rcpp::wrap( kendall_worker_d(SRatings) );
    } catch( std::exception &ex ) {
      forward_exception_to_r( ex );
    } catch(...) { 
      ::Rf_error( "c++ exception (unknown reason)" ); 
    }
    return R_NilValue;
}

SEXP voxelwise_kendall(SEXP Slist_CorMaps, SEXP SSeedMaps, SEXP Sseeds, SEXP Svoxs)
{
    try {
        Rcpp::List list_CorMaps(Slist_CorMaps);
        index_type nsubs = list_CorMaps.size();
        double *voxs = NUMERIC_DATA(Svoxs);
        double *seeds = NUMERIC_DATA(Sseeds);
        Rcpp::NumericVector tmp(Sseeds);
        index_type nseeds = tmp.size();
        
        SEXP SSubMaps;
        arma::mat SubMaps(1,1);
        const double* old_sptr = SubMaps.memptr();
        
        arma::mat SeedMaps(1,1);
        const double* old_smptr = sbm_to_arma_xd(SSeedMaps, SeedMaps);
        index_type nvoxs = static_cast<index_type>(SeedMaps.n_rows);
        if (nsubs != (index_type)SeedMaps.n_cols)
            Rf_error("number of columns in seedmap incorrect");
        
        arma::vec ks(nseeds);
        
        index_type seedi, subi, voxi, seed, vox;
        double* smap;
        index_type add;
        // Loop through each seed voxel
        for (seedi = 0; seedi < nseeds; ++seedi) {
            seed = static_cast<index_type>(seeds[seedi] - 1);
            // Combine seed maps across subjects for given voxel
            for (subi = 0; subi < nsubs; ++subi)
            {
                PROTECT(SSubMaps = VECTOR_ELT(Slist_CorMaps, subi));
                sbm_to_arma_xd(SSubMaps, SubMaps);
                UNPROTECT(1);
                smap = SeedMaps.colptr(subi);
                add = 0;
                for (voxi = 0; voxi < nvoxs; ++voxi) {
                    vox = static_cast<index_type>(voxs[voxi+add] - 1);
                    if (vox==seed) {
                        add = 1;
                        vox = static_cast<index_type>(voxs[voxi+add] - 1);
                    }
                    smap[voxi] = SubMaps(seedi,vox);
                }
            }
            
            // Get kendall's w
            ks(seedi) = kendall_worker_d(SSeedMaps);
        }
        
        // Scale
        double d_nvoxs = static_cast<double>(nvoxs);
        double d_nsubs = static_cast<double>(nsubs);
        ks = (12*ks*(d_nvoxs-1)) / (pow(d_nsubs,2)*(pow(d_nvoxs,3)-d_nvoxs));
        
        smap = NULL;
        free_arma(SubMaps, old_sptr);
        free_arma(SeedMaps, old_smptr);
        
        return Rcpp::wrap( ks );
    } catch( std::exception &ex ) {
      forward_exception_to_r( ex );
    } catch(...) { 
      ::Rf_error( "c++ exception (unknown reason)" ); 
    }
    return R_NilValue;
}

SEXP voxelwise_kendall3(SEXP Slist_CorMaps, SEXP SSeedMaps, SEXP Sseeds, SEXP Svoxs)
{
    try {
        Rcpp::List list_CorMaps(Slist_CorMaps);
        index_type nsubs = list_CorMaps.size();
        double *voxs = NUMERIC_DATA(Svoxs);
        double *seeds = NUMERIC_DATA(Sseeds);
        Rcpp::NumericVector tmp(Sseeds);
        index_type nseeds = tmp.size();
        
        SEXP SSubMaps;
        arma::mat SubMaps(1,1);
        const double* old_sptr = SubMaps.memptr();
        
        arma::mat SeedMaps(1,1);
        const double* old_smptr = sbm_to_arma_xd(SSeedMaps, SeedMaps);
        index_type nvoxs = static_cast<index_type>(SeedMaps.n_rows);
        if (nsubs != (index_type)SeedMaps.n_cols)
            Rf_error("number of columns in seedmap incorrect");
        
        arma::vec ks(nseeds);
        
        index_type seedi, subi, voxi, seed, vox;
        double* smap;
        index_type add;
        // Loop through each seed voxel
        for (seedi = 0; seedi < nseeds; ++seedi) {
            seed = static_cast<index_type>(seeds[seedi] - 1);
            // Combine seed maps across subjects for given voxel
            for (subi = 0; subi < nsubs; ++subi)
            {
                PROTECT(SSubMaps = VECTOR_ELT(Slist_CorMaps, subi));
                sbm_to_arma_xd(SSubMaps, SubMaps);
                UNPROTECT(1);
                smap = SeedMaps.colptr(subi);
                for (voxi = 0; voxi < nvoxs; ++voxi) {
                    smap[voxi] = SubMaps(seedi,voxi);
                }
            }
            
            // Get kendall's w
            ks(seedi) = kendall_worker_d(SSeedMaps);
        }
        
        // Scale
        double d_nvoxs = static_cast<double>(nvoxs);
        double d_nsubs = static_cast<double>(nsubs);
        ks = (12*ks*(d_nvoxs-1)) / (pow(d_nsubs,2)*(pow(d_nvoxs,3)-d_nvoxs));
        
        smap = NULL;
        free_arma(SubMaps, old_sptr);
        free_arma(SeedMaps, old_smptr);
        
        return Rcpp::wrap( ks );
    } catch( std::exception &ex ) {
      forward_exception_to_r( ex );
    } catch(...) { 
      ::Rf_error( "c++ exception (unknown reason)" ); 
    }
    return R_NilValue;
}

SEXP voxelwise_kendall3_regress(SEXP Slist_CorMaps, SEXP SSeedMaps, SEXP SX, SEXP Sseeds, SEXP Svoxs)
{
    try {
        Rcpp::List list_CorMaps(Slist_CorMaps);
        index_type nsubs = list_CorMaps.size();
        double *voxs = NUMERIC_DATA(Svoxs);
        double *seeds = NUMERIC_DATA(Sseeds);
        Rcpp::NumericVector tmp(Sseeds);
        index_type nseeds = tmp.size();
        
        SEXP SSubMaps;
        arma::mat SubMaps(1,1);
        const double* old_sptr = SubMaps.memptr();
        
        arma::mat SeedMaps(1,1);
        const double* old_smptr = sbm_to_arma_xd(SSeedMaps, SeedMaps);
        index_type nvoxs = static_cast<index_type>(SeedMaps.n_rows);
        if (nsubs != (index_type)SeedMaps.n_cols)
            Rf_error("number of rows in seedmap incorrect");
        
        arma::mat y(nsubs, nvoxs);
        arma::mat vox_ones = arma::ones(1, nsubs);
        arma::mat sy(nvoxs, 1);
        arma::mat X(1,1);
        const double* old_sx = sbm_to_arma_xd(SX, X);
        
        arma::vec ks(nseeds);
        
        index_type seedi, subi, voxi, seed, vox;
        double* smap;
        index_type add;
        // Loop through each seed voxel
        for (seedi = 0; seedi < nseeds; ++seedi) {
            seed = static_cast<index_type>(seeds[seedi] - 1);
            sy.zeros();
            // Combine seed maps across subjects for given voxel
            for (subi = 0; subi < nsubs; ++subi)
            {
                PROTECT(SSubMaps = VECTOR_ELT(Slist_CorMaps, subi));
                sbm_to_arma_xd(SSubMaps, SubMaps);
                UNPROTECT(1);
                for (voxi = 0; voxi < nvoxs; ++voxi) {
                    y(subi,voxi) = SubMaps(seedi,voxi);
                    sy(voxi,0) += SubMaps(seedi,voxi);
                }
            }
            
            // fit model y ~ X
            arma::mat coef = arma::solve(X, y);

            // residuals
            SeedMaps = arma::trans(y - X*coef);
            
            // add back column means
            sy = sy/nsubs;
            SeedMaps = SeedMaps + (sy * vox_ones);
            
            // Get kendall's w
            ks(seedi) = kendall_worker_d(SSeedMaps);
        }
        
        // Scale
        double d_nvoxs = static_cast<double>(nvoxs);
        double d_nsubs = static_cast<double>(nsubs);
        ks = (12*ks*(d_nvoxs-1)) / (pow(d_nsubs,2)*(pow(d_nvoxs,3)-d_nvoxs));
        
        smap = NULL;
        free_arma(SubMaps, old_sptr);
        free_arma(SeedMaps, old_smptr);
        free_arma(X, old_sx);
        
        return Rcpp::wrap( ks );
    } catch( std::exception &ex ) {
      forward_exception_to_r( ex );
    } catch(...) { 
      ::Rf_error( "c++ exception (unknown reason)" ); 
    }
    return R_NilValue;
}

SEXP gcor_worker(SEXP ScorMaps, SEXP Sthresh, SEXP SthreshType, 
                 SEXP Sseeds, SEXP Svoxs) 
{
    try {
        arma::mat corMaps(1,1);
        const double* old_cptr = sbm_to_arma_xd(ScorMaps, corMaps);
        index_type nr = static_cast<index_type>(corMaps.n_rows);
        index_type nc = static_cast<index_type>(corMaps.n_cols);

        double thresh = DOUBLE_DATA(Sthresh)[0];
        int threshType = INTEGER_DATA(SthreshType)[0];
        
        double *seeds = NUMERIC_DATA(Sseeds);
        double *voxs = NUMERIC_DATA(Svoxs);

        index_type i, j, n, seed, vox;
        LDOUBLE s = 0;
        double tmp;
        arma::vec res(nr);

        for (i = 0; i < nr; ++i) {
            seed = static_cast<index_type>(seeds[i] - 1);
            s = n = 0;

            for (j = 0; j < nc; ++j) {
                vox = static_cast<index_type>(voxs[j] - 1);
                tmp = corMaps(i,j);
                
                if (seed != vox) {
                    switch (threshType) {
                        case 0:
                            s += static_cast<LDOUBLE>(tmp);
                            n += 1;
                            break;
                        case 1:
                        case 2:
                            if (tmp > thresh) {
                                s += static_cast<LDOUBLE>(tmp);
                                n += 1;
                            }
                            break;
                        case 3:
                        case 4:
                            if (tmp < thresh) {
                                s += static_cast<LDOUBLE>(tmp);
                                n += 1;
                            }
                            break;
                        default:
                            Rf_error("unrecognized threshType");
                            break;
                    }
                }
            }
            
            switch (threshType) {
                case 0:
                case 1:
                case 3:
                    res(i) = static_cast<double>(s/static_cast<LDOUBLE>(n));
                    break;
                case 2:
                case 4:
                    res(i) = static_cast<double>(n);
                    break;
                default:
                    Rf_error("unrecognized threshType");
                    break;            
            }
            
            //printf("i = %i; seed = %i; s = %.3f; n = %i; r = %.3f; val = %.3f\n", 
            //        i, seed, s, n, res(i), tmp);
        }

        free_arma(corMaps, old_cptr);

        return Rcpp::wrap( res );
        
    } catch( std::exception &ex ) {
      forward_exception_to_r( ex );
    } catch(...) { 
      ::Rf_error( "c++ exception (unknown reason)" ); 
    }
    
    return R_NilValue;
}
