#ifndef PRIMESVDSOLVER_H
#define PRIMESVDSOLVER_H

#include "NumericalMethods/SVDSolver.H"
#include "NumericalMethods/SparseSVDSolver.H"

#include "oml/vector.h"
#include "oml/diagonalmatrix.h"

struct primme_svds_params;
//
//  Client code can derive from this and control the matvec multiplication operations
//
/*template <class T> class PrimeSVDSolverClient
{
public:
    virtual int  GetNumRows() const=0;
    virtual int  GetNumCols() const=0;
    virtual void DoMatVecContraction (int N, T* xvec, T* yvec,int transpose) const=0;

    static const PrimeSVDSolverClient* theClient;
};*/
//
//  Factor dense or sparse matrices using the PRIMME library
//  For efficiency client code should keep an instance of PrimeSVDSolver alive
//  between iterations so that initial guesses get used.
//
template <class T> class PrimeSVDSolver
    : public virtual SVDSolver<T>
    , public virtual SparseSVDSolver<T>
{
    using MatrixT         =typename SVDSolver      <T>::MatrixT;
    using DiagonalMatrixRT=typename SVDSolver      <T>::DiagonalMatrixRT;
    using UsVType         =typename SVDSolver      <T>::UsVType;
    using ClientT         =typename SparseSVDSolver<T>::ClientT;
    using SparseMatrixT   =typename SparseSVDSolver<T>::SparseMatrixT;
    using VectorT         =Vector<T>;
    using SparseSVDSolver<T>::theDenseMatrix;
    using SparseSVDSolver<T>::theSparseMatrix;
public:
    PrimeSVDSolver();
    ~PrimeSVDSolver();

    // Get lowest N eigen values and vectors with no initial guess
    virtual UsVType Solve   (const MatrixT&      ,double eps, int NumSingluarValues);
    virtual UsVType SolveAll(const MatrixT&      ,double eps                       );
    virtual UsVType Solve   (const SparseMatrixT&,double eps, int NumSingularValues);
    virtual UsVType Solve   (const ClientT*      ,double eps, int NumSingularValues);

private:
    UsVType Solve(primme_svds_params&);


    DiagonalMatrixRT itsSingularValues;
    VectorT          itsSingularVectors; //Packed storage of U and VT (or Vdagger)
    int              itsNumGuesses;
};

#endif // PRIMESVDSOLVER_H
