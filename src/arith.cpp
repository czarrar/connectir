#include "connectir/connectir.h"

// Return row sum of big matrix (must be of type double)
// Y = a*X + b
SEXP big_add_multiply_scalar(SEXP SX, SEXP SY,  
                             SEXP Sa, SEXP Sb, 
                             SEXP SX_firstCol, SEXP SX_lastCol, 
                             SEXP SY_firstCol, SEXP SY_lastCol) 
{
    try {     
        double a = DOUBLE_DATA(Sa)[0];
        double b = DOUBLE_DATA(Sb)[0];
        
        arma::mat X;
        sub_sbm_to_arma_xd(SX, X, SX_firstCol, SX_lastCol);
        arma::mat Y;
        sub_sbm_to_arma_xd(SY, Y, SY_firstCol, SY_lastCol);
        
        Y = a*X + b;
        
        return R_NilValue;
    } catch(std::exception &ex) {
        forward_exception_to_r(ex);
    } catch(...) {
        ::Rf_error("c++ exception (unknown reason)");
    }    
    return R_NilValue;
}
