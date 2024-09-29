#ifndef SPARSEEIGENSOLVER_H
#define SPARSEEIGENSOLVER_H

#include <tuple>

template <class T> class SparseMatrix;
template <class T> class Matrix;
template <class T> class Vector;
namespace std
{
template <class T> class complex;
}

//
//  Client code can derive from this and control the matvec multiplication operations
//
template <class T> class SparseEigenSolverClient
{
public:
    virtual long int GetSize  () const=0; //2^L Hilbert space size
    virtual void     DoMatVecContraction (int N, const T* xvec, T* yvec) const=0;

    static const SparseEigenSolverClient* theClient;
};

//
//  Client for non-symmetric
//
template <class T> class SparseEigenSolverClientNS
{
public:
    virtual long int GetSize  () const=0; //2^L Hilbert space size
    virtual void     DoMatVecContraction (int N, const T* xvec, T* yvec) const=0;
    virtual void     DoMatDaggerVecContraction (int N, const T* xvec, T* yvec) const=0;

    static const SparseEigenSolverClientNS* theClient;
};
//
//  Abstract interface for eigen solvers
//
template <class T> class SparseEigenSolver
{
protected:
    typedef std::complex<double> dcmplx;
    typedef Vector <double> VectorR;
    typedef Vector <dcmplx> VectorC;
    typedef SparseMatrix<T> SparseMatrixT;
    typedef Matrix<T>       MatrixT;
    typedef Matrix<dcmplx>  MatrixC;
    typedef std::tuple<MatrixT,VectorR> UdType;
    typedef std::tuple<MatrixC,VectorC> UdTypeN;
    typedef SparseEigenSolverClient  <T> ClientT;
    typedef SparseEigenSolverClientNS<T> ClientNST;
public:
    virtual ~SparseEigenSolver() {};
    virtual void Reset()=0;  //Clear out guesses from previous iterations

    virtual UdType Solve(const SparseMatrixT&,double eps, int NumEigenValues)=0;
    virtual UdType Solve(const ClientT*      ,double eps, int NumEigenValues)=0;

    virtual UdTypeN SolveRightNonSym(const SparseMatrixT&,double eps, int NumEigenValues)=0;
    virtual UdTypeN SolveRightNonSym(const ClientT*      ,double eps, int NumEigenValues)=0;
    virtual UdTypeN SolveLeft_NonSym(const SparseMatrixT&,double eps, int NumEigenValues);
//    virtual UdTypeN SolveLeft_NonSym(const ClientT*      ,double eps, int NumEigenValues);

//  No SolveAll here, one purpose of sparse is to not solve all.
    static const SparseMatrix<T>* theSparseMatrix;
    static const       Matrix<T>* theDenseMatrix;

};

#endif // SPARSEEIGENSOLVER_H
