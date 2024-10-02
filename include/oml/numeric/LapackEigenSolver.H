#ifndef LAPACKEIGENSOLVER_H
#define LAPACKEIGENSOLVER_H

#include "oml/numeric/EigenSolver.H"

namespace oml
{

template <class T> class LapackEigenSolver : public virtual EigenSolver<T>
{
    using UdType =typename EigenSolver<T>::UdType;
    using UdTypeN=typename EigenSolver<T>::UdTypeN;
    using MatrixT=typename EigenSolver<T>::MatrixT;
public:
    ~LapackEigenSolver() {};
    virtual void Reset() {};  //Noop Clear out guesses from previous iterations

    virtual UdType  Solve              (const MatrixT&,double eps, int NumEigenValues);
    virtual UdType  SolveAll           (const MatrixT&,double eps                    );
    virtual UdTypeN SolveRightNonSym   (const MatrixT&,double eps, int NumEigenValues);
    virtual UdTypeN SolveAllRightNonSym(const MatrixT&,double eps                    );
};

}
#endif // LAPACKEIGENSOLVER_H
