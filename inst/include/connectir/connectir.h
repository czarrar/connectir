#ifndef _connectir_CONNECTIR
#define _connectir_CONNECTIR

#include <RcppArmadillo.h>
//using namespace Rcpp; 
//using namespace Rcpp::sugar;

#include "connectirDefines.h"

#include <math.h>
#include <iostream>

#include "bigmemory/BigMatrix.h"
#include "bigmemory/MatrixAccessor.hpp"
#include "bigmemory/bigmemoryDefines.h"
#include "bigmemory/isna.hpp"

#include <R.h>
#include <Rinternals.h>
#include <Rdefines.h>

extern "C" {
    SEXP ComputePvalsMain(SEXP Rinmat, SEXP Routmat, SEXP Routcol);
}

RcppExport SEXP CombineSubMapsMain(SEXP LIST_allVoxs_allSubs, SEXP ADDR_oneVox_allSubs, SEXP Rseed_index, SEXP Rvoxindices, SEXP Rnvoxs, SEXP Rnsubs);

RcppExport SEXP CombineSubMapsTransSimpleMain(SEXP LIST_allVoxs_allSubs, SEXP ADDR_oneVox_allSubs, SEXP Rseed_index, SEXP Rvoxindices, SEXP Rnvoxs, SEXP Rnsubs);

// utils.cpp
void free_arma(arma::mat& A, const double *ptr_double);
BigMatrix* sbm_to_bm(SEXP Sbm);
double* bm_to_ptr_xd(BigMatrix* pMat);
const double* sbm_to_arma_xd(SEXP SbM, arma::mat& M);
const double* sub_sbm_to_arma_xd(SEXP SbM, arma::mat& M, SEXP SfirstCol, SEXP SlastCol);
RcppExport SEXP test_func(SEXP Sbm);

// arith.cpp
RcppExport SEXP big_add_multiply_scalar(SEXP SX, SEXP SY,  
                                        SEXP Sa, SEXP Sb, 
                                        SEXP SX_firstCol, SEXP SX_lastCol, 
                                        SEXP SY_firstCol, SEXP SY_lastCol);

// summary_stats.cpp
RcppExport SEXP big_rowsum(SEXP SX, SEXP SX_firstCol, SEXP SX_lastCol);
RcppExport SEXP big_rowmean(SEXP SX, SEXP SX_firstCol, SEXP SX_lastCol);
RcppExport SEXP big_colsum(SEXP SX, SEXP SX_firstCol, SEXP SX_lastCol);
RcppExport SEXP big_colmean(SEXP SX, SEXP SX_firstCol, SEXP SX_lastCol);

// qlm.cpp
RcppExport SEXP big_qlm_rank(SEXP Xr);
RcppExport SEXP big_qlm_dd(SEXP Xr);
RcppExport SEXP big_qlm_fit(SEXP yr, SEXP Xr, SEXP coefr, SEXP residr, SEXP mser, SEXP nr, SEXP kr, SEXP mr);
RcppExport SEXP big_qlm_rsquared(SEXP yr, SEXP coefr, SEXP residr, SEXP mser, SEXP nr, SEXP kr, SEXP rsquaredr, SEXP toAdjust);
RcppExport SEXP big_qlm_residuals(SEXP yr, SEXP Xr, SEXP residr, SEXP add_meanr);
RcppExport SEXP big_qlm_contrasts(SEXP fit_coefr, SEXP fit_mser, SEXP conr, SEXP ddr, SEXP coefr, SEXP ser, SEXP tvalr, SEXP mr);

// cor.cpp
RcppExport SEXP test_sub_matrix(SEXP As, SEXP As_firstCol, SEXP As_lastCol);
RcppExport SEXP big_cor(SEXP As, SEXP Bs, SEXP Cs, 
                        SEXP A_firstCol, SEXP A_lastCol, 
                        SEXP B_firstCol, SEXP B_lastCol, 
                        SEXP C_firstCol, SEXP C_lastCol);
RcppExport SEXP big_tcor(SEXP As, SEXP Bs, SEXP Cs, 
                         SEXP A_firstCol, SEXP A_lastCol, 
                         SEXP B_firstCol, SEXP B_lastCol, 
                         SEXP C_firstCol, SEXP C_lastCol);
RcppExport SEXP kendall_worker(SEXP SRatings);
RcppExport SEXP voxelwise_kendall(SEXP Slist_CorMaps, SEXP SSeedMaps, 
                                  SEXP Sseeds, SEXP Svoxs);
RcppExport SEXP voxelwise_kendall3(SEXP Slist_CorMaps, SEXP SSeedMaps, 
                                   SEXP Sseeds, SEXP Svoxs);
RcppExport SEXP voxelwise_kendall3_regress(SEXP Slist_CorMaps, SEXP SSeedMaps, 
                                           SEXP SX, SEXP Sseeds, SEXP Svoxs);
RcppExport SEXP gcor_worker(SEXP ScorMaps, SEXP Sthresh, SEXP SthreshType, 
                            SEXP Sseeds, SEXP Svoxs);

// subdist2.cpp
RcppExport SEXP subdist_combine_and_scale_submaps(SEXP Slist_corMaps, SEXP Sseed, 
                                                  SEXP SvoxInds, SEXP SseedCorMaps);
RcppExport SEXP subdist_combine_submaps(SEXP Slist_corMaps, SEXP Sseed, 
                                        SEXP SvoxInds, SEXP SseedCorMaps);
RcppExport SEXP subdist_combine_and_trans_submaps(SEXP Slist_corMaps, SEXP Sseed, 
                                                  SEXP SvoxInds, SEXP SseedCorMaps);
RcppExport SEXP subdist_pearson_distance(SEXP SseedCorMaps, SEXP Sdmats, SEXP Sdcol, 
                                         SEXP trans);
RcppExport SEXP big_gower(SEXP SX, SEXP SY,  
                          SEXP SX_firstCol, SEXP SX_lastCol, 
                          SEXP SY_firstCol, SEXP SY_lastCol);

// adonis.cpp
RcppExport SEXP mdmr_fstats_to_pvals(SEXP SFmat);
RcppExport SEXP mdmr_fstats_to_pvals2(SEXP SFperms, SEXP Snperms, SEXP SPmat);
RcppExport SEXP mdmr_worker(SEXP SGmat, SEXP SFperms, 
                            SEXP SH2mats, SEXP SIHmat, 
                            SEXP SdfRes, SEXP SdfExp);


#endif // _connectir_CONNECTIR
