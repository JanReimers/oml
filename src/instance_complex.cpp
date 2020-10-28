//---------------------------------------------------------------------------------
//
//  Make template instance
//
#include "src/matrix.cpp"
#include "src/smatrix.cpp"
#include "src/diagonalmatrix.cpp"
#include "src/vector.cpp"
#define Type std::complex<double>

template class DiagonalMatrix<Type>;
template class         Matrix<Type>;
template class        SMatrix<Type>;
template class         Vector<Type>;

