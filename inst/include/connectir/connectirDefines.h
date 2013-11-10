#ifndef _connectir_DEFINES_H
#define _connectir_DEFINES_H

// take an 'input' R big matrix variable to an 'output' double armadillo matrix
#define BM_TO_ARMA_ONCE(INPUT, OUTPUT) \
    Rcpp::RObject bm(INPUT); \
    SEXP addr = bm.slot("address"); \
    BigMatrix *pMat = reinterpret_cast<BigMatrix*>(R_ExternalPtrAddr(addr)); \
    if (pMat->matrix_type() != 8) \
        ::Rf_error("Big Matrix must be of type double"); \
    index_type offset = pMat->nrow() * pMat->col_offset(); \
    double *ptr_double = reinterpret_cast<double*>(pMat->matrix()) + offset; \
    arma::mat OUTPUT(ptr_double, pMat->nrow(), pMat->ncol(), false);

#define BM_COL_INIT(INPUT, OUTPUT)  \
    index_type OUTPUT = static_cast<index_type>(DOUBLE_DATA(INPUT)[0] - 1);

#define BM_TO_ARMA_INIT()   \
    Rcpp::RObject bm;       \
    SEXP addr;              \
    BigMatrix *pMat;        \
    index_type offset;      \
    index_type ncol;        \
    double *ptr_double;

#define BM_TO_ARMA_MULTIPLE(INPUT, OUTPUT) \
    bm = INPUT; \
    addr = bm.slot("address"); \
    pMat = reinterpret_cast<BigMatrix*>(R_ExternalPtrAddr(addr)); \
    if (pMat->matrix_type() != 8) \
        ::Rf_error("Big Matrix must be of type double"); \
    offset = pMat->nrow() * pMat->col_offset(); \
    ptr_double = reinterpret_cast<double*>(pMat->matrix()) + offset; \
    ncol = pMat->ncol() - pMat->col_offset(); \
    arma::mat OUTPUT(ptr_double, pMat->nrow(), ncol, false);

#define SUB_BM_TO_ARMA_MULTIPLE(INPUT, OUTPUT, COL_FIRST, COL_LAST) \
    if (COL_LAST < COL_FIRST) \
        ::Rf_error("last column smaller than first column"); \
    bm = INPUT; \
    addr = bm.slot("address"); \
    pMat = reinterpret_cast<BigMatrix*>(R_ExternalPtrAddr(addr)); \
    if (pMat->matrix_type() != 8) \
        ::Rf_error("Big Matrix must be of type double"); \
    offset = pMat->nrow() * (pMat->col_offset() + COL_FIRST); \
    ptr_double = reinterpret_cast<double*>(pMat->matrix()) + offset; \
    ncol = COL_LAST - pMat->col_offset() - COL_FIRST + 1; \
    arma::mat OUTPUT(ptr_double, pMat->nrow(), ncol, false);

#define SET_ACCESSOR(ptr, mat)                                          \
    if (ptr->separated_columns()) {                                     \
        switch (ptr->matrix_type()) {                                   \
            case 1:                                                     \
                SepMatrixAccessor<char> mat( *ptr );                    \
                break;                                                  \
            case 2:                                                     \
                SepMatrixAccessor<short> mat( *ptr );                   \
                break;                                                  \
            case 4:                                                     \
                SepMatrixAccessor<int> mat( *ptr );                     \
                break;                                                  \
            case 8:                                                     \
                SepMatrixAccessor<double> mat( *ptr );                  \
                break;                                                  \
        }                                                               \
    }                                                                   \
    else {                                                              \
        switch (ptr->matrix_type()) {                                   \
            case 1:                                                     \
                MatrixAccessor<char> mat( *ptr );                       \
                break;                                                  \
            case 2:                                                     \
                MatrixAccessor<short> mat( *ptr );                      \
                break;                                                  \
            case 4:                                                     \
                MatrixAccessor<int> mat( *ptr );                        \
                break;                                                  \
            case 8:                                                     \
                MatrixAccessor<double> mat( *ptr );                     \
                break;                                                  \
        }                                                               \
    }

#define CALL_BIGFUNCTION_ARGS_ONE(FUN, INBIGMAT)                        \
    SEXP ret;                                                           \
    if (INBIGMAT->separated_columns()) {                                \
        switch (INBIGMAT->matrix_type()) {                              \
            case 1:                                                     \
                ret = FUN<char, SepMatrixAccessor<char> >(                    \
                INBIGMAT);                                              \
                break;                                                  \
            case 2:                                                     \
                ret = FUN<short, SepMatrixAccessor<short> >(                  \
                INBIGMAT);                                              \
                break;                                                  \
            case 4:                                                     \
                ret = FUN<int, SepMatrixAccessor<int> >(                      \
                INBIGMAT);                                              \
                break;                                                  \
            case 8:                                                     \
                ret = FUN<double, SepMatrixAccessor<double> >(                \
                INBIGMAT);                                              \
                break;                                                  \
        }                                                               \
    }                                                                   \
    else {                                                              \
        switch (INBIGMAT->matrix_type()) {                              \
            case 1:                                                     \
                ret = FUN<char, MatrixAccessor<char> >(                       \
                INBIGMAT);                                              \
                break;                                                  \
            case 2:                                                     \
                ret = FUN<short, MatrixAccessor<short> >(                     \
                INBIGMAT);                                              \
                break;                                                  \
            case 4:                                                     \
                ret = FUN<int, MatrixAccessor<int> >(                         \
                INBIGMAT);                                              \
                break;                                                  \
            case 8:                                                     \
                ret = FUN<double, MatrixAccessor<double> >(                   \
                INBIGMAT);                                              \
                break;                                                  \
        }                                                               \
    }

#define CALL_BIGFUNCTION_ARGS_TWO(FUN, INBIGMAT, ARG2)                  \
    SEXP ret;                                                           \
    if (INBIGMAT->separated_columns()) {                                \
        switch (INBIGMAT->matrix_type()) {                              \
            case 1:                                                     \
                ret = FUN<char, SepMatrixAccessor<char> >(                    \
                INBIGMAT, ARG2);                                        \
                break;                                                  \
            case 2:                                                     \
                ret = FUN<short, SepMatrixAccessor<short> >(                  \
                INBIGMAT, ARG2);                                        \
                break;                                                  \
            case 4:                                                     \
                ret = FUN<int, SepMatrixAccessor<int> >(                      \
                INBIGMAT, ARG2);                                        \
                break;                                                  \
            case 8:                                                     \
                ret = FUN<double, SepMatrixAccessor<double> >(                \
                INBIGMAT, ARG2);                                        \
                break;                                                  \
        }                                                               \
    }                                                                   \
    else {                                                              \
        switch (INBIGMAT->matrix_type()) {                              \
            case 1:                                                     \
                ret = FUN<char, MatrixAccessor<char> >(                       \
                INBIGMAT, ARG2);                                        \
                break;                                                  \
            case 2:                                                     \
                ret = FUN<short, MatrixAccessor<short> >(                     \
                INBIGMAT, ARG2);                                        \
                break;                                                  \
            case 4:                                                     \
                ret = FUN<int, MatrixAccessor<int> >(                         \
                INBIGMAT, ARG2);                                        \
                break;                                                  \
            case 8:                                                     \
                ret = FUN<double, MatrixAccessor<double> >(                   \
                INBIGMAT, ARG2);                                        \
                break;                                                  \
        }                                                               \
    }

#define CALL_BIGFUNCTION_ARGS_THREE(FUN, INBIGMAT, ARG2, ARG3)                  \
    SEXP ret;                                                           \
    if (INBIGMAT->separated_columns()) {                                \
        switch (INBIGMAT->matrix_type()) {                              \
            case 1:                                                     \
                ret = FUN<char, SepMatrixAccessor<char> >(                    \
                INBIGMAT, ARG2, ARG3);                                        \
                break;                                                  \
            case 2:                                                     \
                ret = FUN<short, SepMatrixAccessor<short> >(                  \
                INBIGMAT, ARG2, ARG3);                                        \
                break;                                                  \
            case 4:                                                     \
                ret = FUN<int, SepMatrixAccessor<int> >(                      \
                INBIGMAT, ARG2, ARG3);                                        \
                break;                                                  \
            case 8:                                                     \
                ret = FUN<double, SepMatrixAccessor<double> >(                \
                INBIGMAT, ARG2, ARG3);                                        \
                break;                                                  \
        }                                                               \
    }                                                                   \
    else {                                                              \
        switch (INBIGMAT->matrix_type()) {                              \
            case 1:                                                     \
                ret = FUN<char, MatrixAccessor<char> >(                       \
                INBIGMAT, ARG2, ARG3);                                        \
                break;                                                  \
            case 2:                                                     \
                ret = FUN<short, MatrixAccessor<short> >(                     \
                INBIGMAT, ARG2, ARG3);                                        \
                break;                                                  \
            case 4:                                                     \
                ret = FUN<int, MatrixAccessor<int> >(                         \
                INBIGMAT, ARG2, ARG3);                                        \
                break;                                                  \
            case 8:                                                     \
                ret = FUN<double, MatrixAccessor<double> >(                   \
                INBIGMAT, ARG2, ARG3);                                        \
                break;                                                  \
        }                                                               \
    }

#define CALL_BIGFUNCTION_ARGS_FOUR(FUN, INBIGMAT, ARG2, ARG3, ARG4)                  \
    SEXP ret;                                                           \
    if (INBIGMAT->separated_columns()) {                                \
        switch (INBIGMAT->matrix_type()) {                              \
            case 1:                                                     \
                ret = FUN<char, SepMatrixAccessor<char> >(                    \
                INBIGMAT, ARG2, ARG3, ARG4);                                        \
                break;                                                  \
            case 2:                                                     \
                ret = FUN<short, SepMatrixAccessor<short> >(                  \
                INBIGMAT, ARG2, ARG3, ARG4);                                        \
                break;                                                  \
            case 4:                                                     \
                ret = FUN<int, SepMatrixAccessor<int> >(                      \
                INBIGMAT, ARG2, ARG3, ARG4);                                        \
                break;                                                  \
            case 8:                                                     \
                ret = FUN<double, SepMatrixAccessor<double> >(                \
                INBIGMAT, ARG2, ARG3, ARG4);                                        \
                break;                                                  \
        }                                                               \
    }                                                                   \
    else {                                                              \
        switch (INBIGMAT->matrix_type()) {                              \
            case 1:                                                     \
                ret = FUN<char, MatrixAccessor<char> >(                       \
                INBIGMAT, ARG2, ARG3, ARG4);                                        \
                break;                                                  \
            case 2:                                                     \
                ret = FUN<short, MatrixAccessor<short> >(                     \
                INBIGMAT, ARG2, ARG3, ARG4);                                        \
                break;                                                  \
            case 4:                                                     \
                ret = FUN<int, MatrixAccessor<int> >(                         \
                INBIGMAT, ARG2, ARG3, ARG4);                                        \
                break;                                                  \
            case 8:                                                     \
                ret = FUN<double, MatrixAccessor<double> >(                   \
                INBIGMAT, ARG2, ARG3, ARG4);                                        \
                break;                                                  \
        }                                                               \
    }

#define CALL_BIGFUNCTION_ARGS_FIVE(FUN, INBIGMAT, ARG2, ARG3, ARG4, ARG5)                  \
    SEXP ret;                                                           \
    if (INBIGMAT->separated_columns()) {                                \
        switch (INBIGMAT->matrix_type()) {                              \
            case 1:                                                     \
                ret = FUN<char, SepMatrixAccessor<char> >(                    \
                INBIGMAT, ARG2, ARG3, ARG4, ARG5);                                        \
                break;                                                  \
            case 2:                                                     \
                ret = FUN<short, SepMatrixAccessor<short> >(                  \
                INBIGMAT, ARG2, ARG3, ARG4, ARG5);                                        \
                break;                                                  \
            case 4:                                                     \
                ret = FUN<int, SepMatrixAccessor<int> >(                      \
                INBIGMAT, ARG2, ARG3, ARG4, ARG5);                                        \
                break;                                                  \
            case 8:                                                     \
                ret = FUN<double, SepMatrixAccessor<double> >(                \
                INBIGMAT, ARG2, ARG3, ARG4, ARG5);                                        \
                break;                                                  \
        }                                                               \
    }                                                                   \
    else {                                                              \
        switch (INBIGMAT->matrix_type()) {                              \
            case 1:                                                     \
                ret = FUN<char, MatrixAccessor<char> >(                       \
                INBIGMAT, ARG2, ARG3, ARG4, ARG5);                                        \
                break;                                                  \
            case 2:                                                     \
                ret = FUN<short, MatrixAccessor<short> >(                     \
                INBIGMAT, ARG2, ARG3, ARG4, ARG5);                                        \
                break;                                                  \
            case 4:                                                     \
                ret = FUN<int, MatrixAccessor<int> >(                         \
                INBIGMAT, ARG2, ARG3, ARG4, ARG5);                                        \
                break;                                                  \
            case 8:                                                     \
                ret = FUN<double, MatrixAccessor<double> >(                   \
                INBIGMAT, ARG2, ARG3, ARG4, ARG5);                                        \
                break;                                                  \
        }                                                               \
    }

#define CALL_BIGFUNCTION_ARGS_SIX(FUN, INBIGMAT, ARG2, ARG3, ARG4, ARG5, ARG6)                  \
    SEXP ret;                                                           \
    if (INBIGMAT->separated_columns()) {                                \
        switch (INBIGMAT->matrix_type()) {                              \
            case 1:                                                     \
                ret = FUN<char, SepMatrixAccessor<char> >(                    \
                INBIGMAT, ARG2, ARG3, ARG4, ARG5, ARG6);                                        \
                break;                                                  \
            case 2:                                                     \
                ret = FUN<short, SepMatrixAccessor<short> >(                  \
                INBIGMAT, ARG2, ARG3, ARG4, ARG5, ARG6);                                        \
                break;                                                  \
            case 4:                                                     \
                ret = FUN<int, SepMatrixAccessor<int> >(                      \
                INBIGMAT, ARG2, ARG3, ARG4, ARG5, ARG6);                                        \
                break;                                                  \
            case 8:                                                     \
                ret = FUN<double, SepMatrixAccessor<double> >(                \
                INBIGMAT, ARG2, ARG3, ARG4, ARG5, ARG6);                                        \
                break;                                                  \
        }                                                               \
    }                                                                   \
    else {                                                              \
        switch (INBIGMAT->matrix_type()) {                              \
            case 1:                                                     \
                ret = FUN<char, MatrixAccessor<char> >(                       \
                INBIGMAT, ARG2, ARG3, ARG4, ARG5, ARG6);                                        \
                break;                                                  \
            case 2:                                                     \
                ret = FUN<short, MatrixAccessor<short> >(                     \
                INBIGMAT, ARG2, ARG3, ARG4, ARG5, ARG6);                                        \
                break;                                                  \
            case 4:                                                     \
                ret = FUN<int, MatrixAccessor<int> >(                         \
                INBIGMAT, ARG2, ARG3, ARG4, ARG5, ARG6);                                        \
                break;                                                  \
            case 8:                                                     \
                ret = FUN<double, MatrixAccessor<double> >(                   \
                INBIGMAT, ARG2, ARG3, ARG4, ARG5, ARG6);                                        \
                break;                                                  \
        }                                                               \
    }


#endif // _connectir_DEFINES_H
