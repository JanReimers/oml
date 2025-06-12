#ifndef LINEARSOLVER_H
#define LINEARSOLVER_H

template <class T> class Vector;
template <class T> class Matrix;
//
//  Abstract interface for Linear solvers  Ax=b.
//  The combinatorics of Matrix shapes and data types can be overwhelming
//  So we just add things as we need them.
//
template <class T> class LinearSolver
{
protected:
    typedef Matrix<T> MatrixT;
    typedef Vector<T> VectorT;
public:
    virtual ~LinearSolver() {};
//
//  Non packed triangular systems
//      solve A*x=b
//
    virtual VectorT SolveUpperTri(const MatrixT& A,const VectorT& b)=0; // A is upper triangular
    virtual VectorT SolveLowerTri(const MatrixT& A,const VectorT& b)=0; // A is lower triangular
//
//      solve b=x*A
//
//    virtual VectorT SolveUpperTri(const VectorT& b,const MatrixT& A); // A is upper triangular
//    virtual VectorT SolveLowerTri(const VectorT& b,const MatrixT& A); // A is lower triangular
//
//  Full solver
//
    virtual VectorT Solve(const MatrixT& A,const VectorT& b)=0; // A*x=b
//    virtual VectorT Solve(const VectorT& b,const MatrixT& A);   // b=x*A
};

#endif // LINEARSOLVER_H
