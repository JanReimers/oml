//---------------------------------------------------------------------------------
//
//  Make template instance
//
#include "src/dmatrix.cc"
#include "src/vector.cc"

//---------------------------------------------------------------------------------
//
//  Make template instance
//

template class DMatrix<double>;
template class  Vector<double>;

#define Type double
typedef DMatrix<Type> Mat;
const Store MatStore=Full;
#include "oml/imp/matsub.ci"



