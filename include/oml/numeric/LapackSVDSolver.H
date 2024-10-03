#ifndef LAPACKSVDSOLVER_H
#define LAPACKSVDSOLVER_H

#include "oml/numeric/SVDSolver.H"

namespace oml {

template <class T> class LapackSVDSolver
    : public virtual SVDSolver<T>
{
    using MatrixT         =typename SVDSolver<T>::MatrixT;
    using DiagonalMatrixRT=typename SVDSolver<T>::DiagonalMatrixRT;
    using UsVType         =typename SVDSolver<T>::UsVType;
 public:
    LapackSVDSolver() {};
    ~LapackSVDSolver() {};

    virtual UsVType  Solve   (const MatrixT&, int NumSingluarValues);
    virtual UsVType  SolveAll(const MatrixT&                       );
};

}
#endif // LAPACKSVDSOLVER_H
