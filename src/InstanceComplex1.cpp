//---------------------------------------------------------------------------------
//
//  Make template instance
//
#include "oml/src/dmatrix.cc"

typedef std::complex<double> dcmplx;
template class DMatrix<dcmplx>;

#define Type std::complex<double>
typedef DMatrix<Type> Mat;
const Store MatStore=Full;
#include "oml/matsub.ci"



