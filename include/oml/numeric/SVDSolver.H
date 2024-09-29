#ifndef SVDSOLVER_H
#define SVDSOLVER_H

#include <tuple>

template <class T> class Matrix;
template <class T> class DiagonalMatrix;
//
//  Abstract interface for SVD solvers
//
template <class T> class SVDSolver
{
protected:
    typedef         Matrix<T>               MatrixT;
    typedef DiagonalMatrix<double> DiagonalMatrixRT;
    typedef std::tuple<MatrixT,DiagonalMatrixRT,MatrixT> UsVType;
public:
    virtual ~SVDSolver() {};

    virtual UsVType  Solve   (const MatrixT&,double eps, int NumSingluarValues)=0;
    virtual UsVType  SolveAll(const MatrixT&,double eps                       )=0;

};


#include "TensorNetworks/Enums.H"
#include <vector>
#include "oml/imp/index_t.h"
//
//  In order to sort the singular values in descending order (required for compression)canned SVD
//  routines must shuffle the rows and columns of U and VT.  THis destroys triangularity.
//  These routines attempt to find a re-indexing (row or column re-shuffling) the restores
//  triangularity.
//  Create a re-indexing array to bring the matrix UL into Upper/Lower triangular form
//
std::vector<index_t>  FindRowReIndex(TensorNetworks::TriType ul,const Matrix<double>& UL, double eps);
std::vector<index_t>  FindColReIndex(TensorNetworks::TriType ul,const Matrix<double>& UL, double eps);

#endif // SVDSOLVER_H
