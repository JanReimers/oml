#ifndef ARPACKEIGENSOLVER_H
#define ARPACKEIGENSOLVER_H

#include "EigenSolver.H"
#include "SparseEigenSolver.H"

template <class T> class ArpackEigenSolver
  : public virtual SparseEigenSolver<T>
  , public virtual       EigenSolver<T>
{
    using MatrixT=typename EigenSolver<T>::MatrixT;
    using SparseMatrixT=typename SparseEigenSolver<T>::SparseMatrixT;

    using VectorC=typename SparseEigenSolver<T>::VectorC;
    using MatrixC=typename SparseEigenSolver<T>::MatrixC;
    using dcmplx =typename SparseEigenSolver<T>::dcmplx;
    using UdType =typename SparseEigenSolver<T>::UdType;
    using UdTypeN=typename EigenSolver<T>::UdTypeN;
public:
    using ClientT=typename SparseEigenSolver<T>::ClientT;

    virtual ~ArpackEigenSolver() {};
    virtual void Reset() {};  //Clear out guesses from previous iterations

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

    //Return value is: (EigenValues,EigenVectors)
//    UdType SolveNonSym(const Matrix     <T>&, int NumEigenValues,double eps);
//    UdType SolveNonSym(const SparseMatrix<T>&, int NumEigenValues,double eps);
    static const SparseMatrix<T>* theSparseMatrix;
    static const       Matrix<T>* theDenseMatrix;
    typedef    void (*MatvecT) (int N, const T * x, T * y);
private:
    // Generic version
    UdType  SolveSym(MatvecT,int N, int NumEigenValues,double eps);
    UdTypeN SolveG  (MatvecT,int N, int NumEigenValues,double eps);
};

#endif // ARPACKEIGENSOLVER_H
