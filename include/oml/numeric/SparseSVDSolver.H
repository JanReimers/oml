#ifndef SPARSESVDSOLVER_H
#define SPARSESVDSOLVER_H

#include <tuple>

template <class T> class SparseMatrix;
template <class T> class DiagonalMatrix;
template <class T> class Matrix;
template <class T> class Vector;
namespace std
{
template <class T> class complex;
}

//
//  Client code can derive from this and control the matvec multiplication operations
//
template <class T> class SparseSVDSolverClient
{
public:
    virtual int  GetNumRows() const=0;
    virtual int  GetNumCols() const=0;
//    virtual int  GetSize  () const=0; //2^L Hilbert space size
    virtual void DoMatVecContraction (int N, const T* xvec, T* yvec,int transpose) const=0;

    static const SparseSVDSolverClient<T>* theClient;
};
//
//  Abstract interface for eigen solvers
//
template <class T> class SparseSVDSolver
{
protected:
    typedef         Vector<T>                VectorT;
    typedef         Matrix<T>                MatrixT;
    typedef   SparseMatrix<T>          SparseMatrixT;
    typedef DiagonalMatrix<double>  DiagonalMatrixRT;
    typedef SparseSVDSolverClient<T>         ClientT;
    typedef std::tuple<MatrixT,DiagonalMatrixRT,MatrixT> UsVType;
public:
    virtual ~SparseSVDSolver() {};

    virtual UsVType Solve(const SparseMatrixT&,double eps, int NumSingularValues)=0;
    virtual UsVType Solve(const ClientT*      ,double eps, int NumSingularValues)=0;

//  No SolveAll here, only purpose of sparse is to not solve all.
    static const SparseMatrix<T>* theSparseMatrix;
    static const       Matrix<T>* theDenseMatrix;

};

#endif // SPARSESVDSOLVER_H
