//---------------------------------------------------------------------------------
//
//  Make template instance
//
#include "oml/src/dmatrix.cc"
#include "oml/src/vector.cc"

//---------------------------------------------------------------------------------
//
//  Make template instance
//
typedef std::complex<double> dcmplx;

template class DMatrix<double>;
template class  Vector<double>;
template class  Vector<dcmplx>;

#define Type double
typedef DMatrix<Type> Mat;
const Store MatStore=Full;
#include "oml/matsub.ci"



