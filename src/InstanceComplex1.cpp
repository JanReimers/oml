//---------------------------------------------------------------------------------
//
//  Make template instance
//
#include "src/dmatrix.cc"
#include "src/vector.cc"

typedef std::complex<double> dcmplx;
template class DMatrix<dcmplx>;
template class  Vector<dcmplx>;

#define Type std::complex<double>
typedef DMatrix<Type> Mat;
const Store MatStore=Full;
#include "oml/imp/matsub.ci"
