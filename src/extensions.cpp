//#define BIG_FUNCTION_SETUP(invar, outvar)    \
//  BMAccessorType outvar( *invar );     \
//    \
//  index_type i=0;   \
//  index_type j=0;   \
//  index_type numCols = pMat->ncol();    \
//  index_type numRows = pMat->nrow();
//
//// Can one restrict CType to float or double?
//template<typename CType>
//void big_pow(CType *invar, CType exponent, index_type nelem, CType *outvar) {
//    for (index_type i = 0; i < nelem; ++i)
//        outvar[i] = pow(invar, exponent);
//    return;
//}
//
//template<typename CType>
//void big_sqrt(CType *invar, index_type nelem, CType *outvar) {
//    for (index_type i = 0; i < nelem; ++i)
//        outvar[i] = sqrt(invar);
//    return;
//}
//
//template<typename CType>
//void big_atanh(CType *invar, index_type nelem, CType *outvar) {
//    for (index_type i = 0; i < nelem; ++i)
//        outvar[i] = atanh(invar);
//    return;
//}
//
//template<typename CType>
//void big_tanh(CType *invar, index_type nelem, CType *outvar) {
//    for (index_type i = 0; i < nelem; ++i)
//        outvar[i] = tanh(invar);
//    return;
//}
//