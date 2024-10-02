#ifndef EIGENSOLVER_H
#define EIGENSOLVER_H

#include <tuple>

template <class T> class Matrix;
template <class T> class Vector;
namespace std
{
template <class T> class complex;
}

namespace oml {
//
//  Abstract interface for eigen solvers
//
template <class T> class EigenSolver
{
protected:
    typedef std::complex<double> dcmplx;
    typedef Vector <double> VectorR;
    typedef Vector <dcmplx> VectorC;
    typedef Matrix<T>       MatrixT;
    typedef Matrix<dcmplx>  MatrixC;
    typedef std::tuple<MatrixT,VectorR> UdType;
    typedef std::tuple<MatrixC,VectorC> UdTypeN;
public:
    virtual ~EigenSolver() {};
    virtual void Reset()=0;  //Clear out guesses from previous iterations

    virtual UdType  Solve              (const MatrixT&,double eps, int NumEigenValues)=0;
    virtual UdType  SolveAll           (const MatrixT&,double eps                    )=0;
    virtual UdTypeN SolveRightNonSym   (const MatrixT&,double eps, int NumEigenValues)=0;
    virtual UdTypeN SolveLeft_NonSym   (const MatrixT&,double eps, int NumEigenValues)  ;
    virtual UdTypeN SolveAllRightNonSym(const MatrixT&,double eps                    )=0;

};

}

#endif // EIGENSOLVER_H
