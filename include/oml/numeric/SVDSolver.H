#ifndef SVDSOLVER_H
#define SVDSOLVER_H

#include <tuple>

template <class T> class Matrix;
template <class T> class DiagonalMatrix;

namespace oml {
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

    virtual UsVType  Solve   (const MatrixT&, size_t NumSingluarValues)=0;
    virtual UsVType  SolveAll(const MatrixT&                       )=0;

};

}


#endif // SVDSOLVER_H
