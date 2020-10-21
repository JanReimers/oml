//---------------------------------------------------------------------------------
//
//  Make template instances
//
#include "src/matrix.cpp"
#include "src/diagonalmatrix.cpp"
#include "src/vector.cpp"
#define Type double

template class DiagonalMatrix<Type>;
template class         Matrix<Type>;
template class         Vector<Type>;
