#include "PrimeEigenSolver.H"
#include "Containers/SparseMatrix.H"
#include "oml/matrix.h"
#include <primme.h>

using std::cout;
using std::endl;

template <class T, class TM> void  DenseMatvec(void *x, PRIMME_INT *ldx, void *y, PRIMME_INT *ldy, int *blockSize, primme_params *primme, int *ierr);
template <class T, class TM> void SparseMatvec(void *x, PRIMME_INT *ldx, void *y, PRIMME_INT *ldy, int *blockSize, primme_params *primme, int *ierr);
template <class T> void ClientMatvec(void *x, PRIMME_INT *ldx, void *y, PRIMME_INT *ldy, int *blockSize, primme_params *primme, int *ierr);


template <class T> PrimeEigenSolver<T>::PrimeEigenSolver()
: itsNumGuesses(0)
{}

template <class T> PrimeEigenSolver<T>::~PrimeEigenSolver()
{}

template <class T> Vector<T>  PrimeEigenSolver<T>::GetEigenVector(int index) const
{
    return itsEigenVectors.GetColumn(index);
}

// Function pointer type for mat*vec functions
typedef    void (*MatvecT) (void *x, PRIMME_INT *ldx, void *y, PRIMME_INT *ldy, int *blockSize,struct primme_params *primme, int *ierr);
// Function dec for building primme paramater struct
primme_params MakeParameters(MatvecT MatVec,int N,int NumEigenValues,int NumGuesses,double eps);
primme_params MakeParametersNonSym(MatvecT MatVec,int N,int NumEigenValues,int NumGuesses,double eps);

template <class T> typename PrimeEigenSolver<T>::UdType
PrimeEigenSolver<T>::Solve(const MatrixT& m,double eps, int NumEigenValues)
{
    int N=m.GetNumRows();
    assert(m.GetLimits().Row.Low==1);
    assert(m.GetLimits().Col.Low==1);
    assert(N==m.GetNumCols());
    assert(NumEigenValues<=N);
    assert(IsHermitian(m,eps));
    theDenseMatrix=&m;
    primme_params primme=MakeParameters(DenseMatvec<T,T>,N,NumEigenValues,itsNumGuesses,eps);
    return Solve(primme);
}

template <class T> typename PrimeEigenSolver<T>::UdType
PrimeEigenSolver<T>::SolveAll(const MatrixT& A,double eps)
{
    int mn=Min(A.GetNumRows(),A.GetNumCols());
    return Solve(A,eps,mn);
}

template <class T> typename PrimeEigenSolver<T>::UdType
PrimeEigenSolver<T>::Solve(const SparseMatrixT& m,double eps, int NumEigenValues)
{
    int N=m.GetNumRows();
    assert(m.GetLimits().Row.Low==1);
    assert(m.GetLimits().Col.Low==1);
    assert(N==m.GetNumCols());
    assert(NumEigenValues<=N);
    assert(IsHermitian(m,eps));
    theSparseMatrix=&m;
    primme_params primme=MakeParameters(SparseMatvec<T,T>,N,NumEigenValues,itsNumGuesses,eps);
    return Solve(primme);
}

template <class T> typename PrimeEigenSolver<T>::UdType
PrimeEigenSolver<T>::Solve(const ClientT* client,double eps, int NumEigenValues)
{
    assert(client);
    int N=client->GetSize();
    assert(NumEigenValues<=N);
    ClientT::theClient=client;  //Used by the M*v call back
    primme_params primme=MakeParameters(ClientMatvec<T>,N,NumEigenValues,itsNumGuesses,eps);
    return Solve(primme);
}



template <class T> typename PrimeEigenSolver<T>::UdTypeN
PrimeEigenSolver<T>::SolveRightNonSym   (const MatrixT& m,double eps, int NumEigenValues)
{
    int N=m.GetNumRows();
    assert(m.GetLimits().Row.Low==1);
    assert(m.GetLimits().Col.Low==1);
    assert(N==m.GetNumCols());
    assert(NumEigenValues<=N);
    //assert(IsNormal(m,eps));
    theDenseMatrix=&m;
    primme_params primme=MakeParametersNonSym(DenseMatvec<dcmplx,T>,N,NumEigenValues,itsNumGuesses,eps);
    return SolveNormal(primme);
}

template <class T> typename PrimeEigenSolver<T>::UdTypeN
PrimeEigenSolver<T>::SolveAllRightNonSym(const MatrixT& A,double eps)
{
    int N=A.GetNumRows();
    return SolveRightNonSym(A,eps,N);
}

template <class T> typename PrimeEigenSolver<T>::UdTypeN
PrimeEigenSolver<T>::SolveRightNonSym(const SparseMatrixT& m,double eps, int NumEigenValues)
{
    int N=m.GetNumRows();
    assert(m.GetLimits().Row.Low==1);
    assert(m.GetLimits().Col.Low==1);
    assert(N==m.GetNumCols());
    assert(NumEigenValues<=N);
//    assert(IsHermitian(m,eps));
    theSparseMatrix=&m;
    primme_params primme=MakeParametersNonSym(SparseMatvec<dcmplx,T>,N,NumEigenValues,itsNumGuesses,eps);
    return SolveNormal(primme);
}

template <class T> typename PrimeEigenSolver<T>::UdTypeN
PrimeEigenSolver<T>::SolveRightNonSym(const ClientT*      ,double eps, int NumEigenValues)
{
    std::cerr << "PrimeEigenSolver::SolveNonSym for Client is not implemented yet." << std::endl;
    return std::make_tuple(MatrixC(),VectorC());
}

//
//  Enable template solve function to call the correct primme routine
//
typedef std::complex<double> dcmplx;

template <class T> int primmeT(double *evals, T *evecs, double *resNorms, primme_params *primme);
template <> int primmeT<double> (double *evals, double *evecs, double *resNorms, primme_params *primme)
{
    return dprimme(evals,evecs,resNorms,primme); //double
}
template <> int primmeT<dcmplx> (double *evals, dcmplx *evecs, double *resNorms, primme_params *primme)
{
    return zprimme(evals,evecs,resNorms,primme); //complex<double>
}

int zprimme_normal (dcmplx *evals, dcmplx *evecs, double *resNorms, primme_params *primme);

//
//  Lowest level solve routines used by all higher level solve functions
//
template <class T> typename PrimeEigenSolver<T>::UdType
PrimeEigenSolver<T>::Solve(primme_params& p)
{
    itsEigenValues.SetLimits(p.numEvals);
    itsEigenVectors.SetLimits(p.n,p.numEvals);
    Vector<double> rnorms(p.numEvals);
    int ret = primmeT<T>(&itsEigenValues(1), &itsEigenVectors(1,1), &rnorms(1), &p);
    if (ret!=0)
        std::cerr << "Error in primme solver, ret=" << ret << std::endl;
    assert(ret==0);
    (void)ret; //avoid compiler warning in release modems)) << " " << std::endl;
    if (Max(fabs(rnorms))>1000*p.eps)
        cout << "Warning high rnorms in PrimeEigenSolver::SolveSparse Max(rnorms)=" << std::scientific << Max(rnorms) << endl;
    //int niter=p.stats.numOuterIterations;
    //std::cout << "Primme niter=" << niter << std::endl;
    itsNumGuesses=p.numEvals; //Set up using guesses for next time around
    primme_free(&p);
    return std::make_tuple(itsEigenVectors,itsEigenValues);
}

template <class T> typename PrimeEigenSolver<T>::UdTypeN
PrimeEigenSolver<T>::SolveNormal(primme_params& p)
{
    itsNonSymEigenValues.SetLimits(p.numEvals);
    itsNonSymEigenVectors.SetLimits(p.n,p.numEvals);
    Vector<double> rnorms(p.numEvals);
    int ret = zprimme_normal(&itsNonSymEigenValues(1), &itsNonSymEigenVectors(1,1), &rnorms(1), &p);
    if (ret!=0)
        std::cerr << "Error in primme solver, ret=" << ret << std::endl;
    assert(ret==0);
    (void)ret; //avoid compiler warning in release modems)) << " " << std::endl;
    if (Max(fabs(rnorms))>1000*p.eps)
        cout << "Warning high rnorms in PrimeEigenSolver::SolveSparse Max(rnorms)=" << std::scientific << Max(rnorms) << endl;
    //int niter=p.stats.numOuterIterations;
    //std::cout << "Primme niter=" << niter << std::endl;
    itsNumGuesses=p.numEvals; //Set up using guesses for next time around
    primme_free(&p);
    return std::make_tuple(itsNonSymEigenVectors,itsNonSymEigenValues);
}

//
//  Build up parameters structure
//

primme_params MakeParameters(MatvecT MatVec,int N,int NumEigenValues,int NumGuesses,double eps)
{
    primme_params primme;
    primme_initialize(&primme);
    primme.matrixMatvec = MatVec;
    primme.n = N; /* set problem dimension */
    primme.numEvals = NumEigenValues;   /* Number of wanted eigenpairs */
    primme.eps = eps;      /* ||r|| <= eps * ||matrix|| */
    primme.target = primme_smallest; /* Wanted the smallest eigenvalues */

    primme.initSize=NumGuesses;
    primme_set_method(PRIMME_DYNAMIC, &primme);
    return primme;
}
double targetShifts[1];
primme_params MakeParametersNonSym(MatvecT MatVec,int N,int NumEigenValues,int NumGuesses,double eps)
{
    primme_params primme;
    primme_initialize(&primme);
    primme.matrixMatvec = MatVec;
    primme.n = N; /* set problem dimension */
    primme.numEvals = NumEigenValues;   /* Number of wanted eigenpairs */
    primme.eps = eps;      /* ||r|| <= eps * ||matrix|| */
    targetShifts[0] = 1000.0;
    primme.targetShifts = targetShifts;
    primme.numTargetShifts = 1;
    primme.target = primme_closest_abs;//primme_smallest; /* Wanted the smallest eigenvalues */

    primme.initSize=NumGuesses;
    primme.printLevel=1;
    primme_set_method(PRIMME_DEFAULT_MIN_MATVECS, &primme);
    primme.correctionParams.projectors.RightX = 0;
    primme.maxOuterIterations=10000;

    primme_set_method(PRIMME_DYNAMIC, &primme);
    return primme;
}

//  matrix-vector product, Y = A * X
//   - X, input dense matrix of size primme.n x blockSize;
//   - Y, output dense matrix of size primme.n x blockSize;
//   - A, square matrix of dimension primme.n with this form:

template <class T, class TM> void SparseMatvec(void *x, PRIMME_INT *_ldx, void *y, PRIMME_INT *_ldy, int *_blockSize, primme_params *primme, int *ierr)
{
    typedef SparseEigenSolver<TM> SES;
    assert(SES::theSparseMatrix);
    long int& ldx(*_ldx);
    long int& ldy(*_ldy);
    int& blockSize(*_blockSize);

    for (int ib=0; ib<blockSize; ib++)
    {
        T* xvec = static_cast<T*>(x) + ldx*ib;
        T* yvec = static_cast<T*>(y) + ldy*ib;
        SES::theSparseMatrix->DoMVMultiplication(primme->n,xvec,yvec);
    }
    *ierr = 0;
}

template <class T,class TM> void DenseMatvec(void *x, PRIMME_INT *_ldx, void *y, PRIMME_INT *_ldy, int *_blockSize, primme_params *primme, int *ierr)
{
    typedef SparseEigenSolver<TM> SES;
    assert(SES::theDenseMatrix);
    long int& ldx(*_ldx);
    long int& ldy(*_ldy);
    int& blockSize(*_blockSize);
    int N=primme->n;
    for (int ib=0; ib<blockSize; ib++)
    {
        T* xvec = static_cast<T*>(x) + ldx*ib;
        T* yvec = static_cast<T*>(y) + ldy*ib;
        for (int ir=1;ir<=N;ir++)
        {
            yvec[ir-1]=0.0;
            for (int ic=1;ic<=N;ic++)
                yvec[ir-1]+=(*SES::theDenseMatrix)(ir,ic)*xvec[ic-1];
        }
    }
    *ierr = 0;

}



template <class T>
void ClientMatvec(void *x, PRIMME_INT *_ldx, void *y, PRIMME_INT *_ldy, int *_blockSize, primme_params *primme, int *ierr)
{
    typedef typename PrimeEigenSolver<T>::ClientT ClientT;
    assert(ClientT::theClient); //|Psi>
    long int& ldx(*_ldx);
    long int& ldy(*_ldy);
    int& blockSize(*_blockSize);

    for (int ib=0; ib<blockSize; ib++)
    {
        T* xvec = static_cast<T*>(x) + ldx*ib;
        T* yvec = static_cast<T*>(y) + ldy*ib;
        ClientT::theClient->DoMatVecContraction(primme->n,xvec,yvec);
    }
    *ierr = 0;

}


//
//  Make template instances
//
template class PrimeEigenSolver<std::complex<double> >;
template class PrimeEigenSolver<double>;

