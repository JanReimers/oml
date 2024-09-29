#ifndef PRIMEEIGENSOLVER_H
#define PRIMEEIGENSOLVER_H

#include "EigenSolver.H"
#include "SparseEigenSolver.H"
#include "oml/vector.h"

struct primme_params;

//
//  Solve sparse eigen systems using the PRIMME library
//  For efficiency client code should keep and instance of PrimeEigenSolver alive
//  between iterations so that eigen vector initial guess get used.
//
template <class T> class PrimeEigenSolver
  : public virtual SparseEigenSolver<T>
  , public virtual       EigenSolver<T>
{
    using dcmplx =typename EigenSolver<T>::dcmplx;
    using UdType =typename SparseEigenSolver<T>::UdType;
    using UdTypeN=typename EigenSolver<T>::UdTypeN;
    using MatrixT=typename EigenSolver<T>::MatrixT;
    using SparseMatrixT=typename SparseEigenSolver<T>::SparseMatrixT;

    using VectorC=typename SparseEigenSolver<T>::VectorC;
    using MatrixC=typename SparseEigenSolver<T>::MatrixC;

    using SparseEigenSolver<T>::theDenseMatrix;
    using SparseEigenSolver<T>::theSparseMatrix;
public:
    using ClientT=typename SparseEigenSolver<T>::ClientT;

    PrimeEigenSolver();
    ~PrimeEigenSolver();
    virtual void Reset() {itsNumGuesses=0;};  //Clear out guesses from previous iterations
    //Dense interface
    virtual UdType  Solve              (const MatrixT&,double eps, int NumEigenValues);
    virtual UdType  SolveAll           (const MatrixT&,double eps                    );
    virtual UdTypeN SolveRightNonSym   (const MatrixT&,double eps, int NumEigenValues);
    virtual UdTypeN SolveAllRightNonSym(const MatrixT&,double eps                    );
    // Sparse interface
    virtual UdType  Solve      (const SparseMatrixT&,double eps, int NumEigenValues);
    virtual UdType  Solve      (const ClientT*      ,double eps, int NumEigenValues);
    virtual UdTypeN SolveRightNonSym(const SparseMatrixT&,double eps, int NumEigenValues);
    virtual UdTypeN SolveRightNonSym(const ClientT*      ,double eps, int NumEigenValues);

    Vector <T> GetEigenVector (int index) const;
    Matrix<T> GetEigenVectors() const
    {
        return itsEigenVectors;
    }
    Vector <double> GetEigenValues () const
    {
        return itsEigenValues;
    }

protected:
    UdType  Solve(primme_params&);
    UdTypeN SolveNormal(primme_params&);

private:
    Vector<double>   itsEigenValues;
    Matrix<T>        itsEigenVectors;
    Vector<dcmplx>   itsNonSymEigenValues;
    Matrix<dcmplx>   itsNonSymEigenVectors;

    int              itsNumGuesses;
};

#endif // PRIMEEIGENSOLVER_H
