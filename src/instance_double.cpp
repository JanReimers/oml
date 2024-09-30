//---------------------------------------------------------------------------------
//
//  Make template instances
//
#include "src/matrix.cpp"
#include "src/smatrix.cpp"
#include "src/diagonalmatrix.cpp"
#include "src/vector.cpp"
#define Type double

template class DiagonalMatrix<Type>; 
template class         Matrix<Type>;
template class        SMatrix<Type>;
template class         Vector<Type>;

#include "oml/vector3d.h"
//#include "oml/matrix3d.h"
#include "oml/io3d.h"

//template class         Matrix<Matrix3D<Type> >;
template class         Vector<Vector3D<Type> >;
template class         Matrix<Vector3D<Type> >;

