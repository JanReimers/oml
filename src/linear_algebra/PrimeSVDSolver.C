#include "PrimeSVDSolver.H"
#include "Containers/SparseMatrix.H"
#include "oml/matrix.h"
#include <primme.h>

using std::cout;
using std::endl;


template <class T> void  DenseMatvec(void *x, PRIMME_INT *ldx, void *y, PRIMME_INT *ldy, int *blockSize, int* transpose, primme_svds_params *primme, int *ierr);
template <class T> void SparseMatvec(void *x, PRIMME_INT *ldx, void *y, PRIMME_INT *ldy, int *blockSize, int* transpose, primme_svds_params *primme, int *ierr);
template <class T> void ClientMatvec(void *x, PRIMME_INT *ldx, void *y, PRIMME_INT *ldy, int *blockSize, int* transpose, primme_svds_params *primme, int *ierr);


template <class T> PrimeSVDSolver<T>::PrimeSVDSolver()
: itsSingularValues()
, itsSingularVectors(0)
, itsNumGuesses(0)
{}

template <class T> PrimeSVDSolver<T>::~PrimeSVDSolver()
{}

// Function pointer type for mat*vec functions
typedef    void (*MatvecT) (void *x, PRIMME_INT *ldx, void *y, PRIMME_INT *ldy, int *blockSize,int* transpose, struct primme_svds_params *primme, int *ierr);
// Function dec for building primme paramater struct
primme_svds_params MakeParameters(MatvecT MatVec,int M,int N,int NumSingularValues,int NumGuesses,double eps);



template <class T> typename PrimeSVDSolver<T>::UsVType
PrimeSVDSolver<T>::Solve(const MatrixT& m,double eps, int NumSingularValues)
{
    assert(m.GetLimits().Row.Low==1);
    assert(m.GetLimits().Col.Low==1);
    assert(NumSingularValues<=Min(m.GetNumRows(),m.GetNumCols()));
    theDenseMatrix=&m;
    primme_svds_params primme=MakeParameters(DenseMatvec<T>,m.GetNumRows(),m.GetNumCols(),NumSingularValues,itsNumGuesses,eps);
    return Solve(primme);
}

template <class T> typename PrimeSVDSolver<T>::UsVType
PrimeSVDSolver<T>::SolveAll(const MatrixT& A,double eps)
{
    int N=Min(A.GetNumRows(),A.GetNumCols());
    return Solve(A,eps,N);
}

template <class T> typename PrimeSVDSolver<T>::UsVType
PrimeSVDSolver<T>::Solve(const SparseMatrixT& m,double eps, int NumSingularValues)
{
    assert(m.GetLimits().Row.Low==1);
    assert(m.GetLimits().Col.Low==1);
    assert(NumSingularValues<=Min(m.GetNumRows(),m.GetNumCols()));
    theSparseMatrix=&m;
    primme_svds_params primme=MakeParameters(SparseMatvec<T>,m.GetNumRows(),m.GetNumCols(),NumSingularValues,itsNumGuesses,eps);
    return Solve(primme);
}

template <class T> typename PrimeSVDSolver<T>::UsVType
PrimeSVDSolver<T>::Solve(const ClientT* client,double eps, int NumSingularValues)
{
    assert(client);
    assert(NumSingularValues<=Min(client->GetNumRows(),client->GetNumCols()));
    ClientT::theClient=client;  //Used by the M*v call back

    primme_svds_params primme=MakeParameters(ClientMatvec<T>,client->GetNumRows(),client->GetNumCols(),NumSingularValues,itsNumGuesses,eps);
    return Solve(primme);
}

//
//  Enable template solve function to call the correct primme routine
//
typedef std::complex<double> dcmplx;
template <class T> int primmeT(double *svals, T *svecs, double *resNorms, primme_svds_params *primme);
template <> int primmeT<double> (double *svals, double *svecs, double *resNorms, primme_svds_params *primme)
{
    return dprimme_svds(svals,svecs,resNorms,primme); //double
}
template <> int primmeT<dcmplx> (double *svals, dcmplx *svecs, double *resNorms, primme_svds_params *primme)
{
    return zprimme_svds(svals,svecs,resNorms,primme); //complex<double>
}

//
//  Lowest level solve routine used by all higher level solve functions
//
template <class T> typename PrimeSVDSolver<T>::UsVType
PrimeSVDSolver<T>::Solve(primme_svds_params& p)
{
    itsSingularValues.SetLimits(p.numSvals,p.numSvals);
    itsSingularVectors.SetLimits((p.m+p.n)*p.numSvals);
    Vector<double> rnorms(p.numSvals);
    int ret = primmeT<T>(&itsSingularValues(1), &itsSingularVectors(1), &rnorms(1), &p);
    assert(ret==0);
    (void)ret; //avoid compiler warning in release modems)) << " " << std::endl;
    if (Max(fabs(rnorms))>1000*p.eps)
        cout << "Warning high rnorms in PrimeSVDSolver::Solve rnorma=" << std::scientific << rnorms << endl;
//    int niter=p.stats.numOuterIterations;
    //std::cout << "Primme niter=" << niter << std::endl;
    itsNumGuesses=p.numSvals; //Set up using guesses for next time around
    primme_svds_free(&p);
    //
    //  Unpack U and VT
    //

    MatrixT U(p.m,p.numSvals),VT(p.numSvals,p.n);
    int linearIndex=1;
    for (int j=1;j<=p.numSvals;j++)
        for (int i=1;i<=p.m;i++,linearIndex++)
            U(i,j)=itsSingularVectors(linearIndex);
    for (int j=1;j<=p.numSvals;j++)
        for (int i=1;i<=p.n;i++,linearIndex++)
            VT(j,i)=conj(itsSingularVectors(linearIndex));

    return std::make_tuple(std::move(U),std::move(itsSingularValues),std::move(VT));
}

//
//  Build up parameters structure
//
primme_svds_params MakeParameters(MatvecT MatVec,int M,int N,int NumSingularValues,int NumGuesses,double eps)
{
    primme_svds_params primme;
    primme_svds_initialize(&primme);
    primme.matrixMatvec = MatVec;
    primme.m = M; /* set problem dimension */
    primme.n = N; /* set problem dimension */
    primme.numSvals = NumSingularValues;   /* Number of wanted eigenpairs */
    primme.eps = eps;      /* ||r|| <= eps * ||matrix|| */
    primme.target = primme_svds_largest; /* Wanted the smallest eigenvalues */

    primme.initSize=NumGuesses;
    primme_svds_set_method (primme_svds_hybrid, PRIMME_DEFAULT_METHOD,
                         PRIMME_DEFAULT_METHOD, &primme);
    return primme;
}


//  matrix-vector product, Y = A * X
//   - X, input dense matrix of size primme.n x blockSize;
//   - Y, output dense matrix of size primme.n x blockSize;
//   - A, square matrix of dimension primme.n with this form:

template <class T> void SparseMatvec(void *x, PRIMME_INT *_ldx, void *y, PRIMME_INT *_ldy, int *_blockSize,int* transpose,  primme_svds_params *primme, int *ierr)
{
    typedef SparseSVDSolver<T> SSS;
    assert(SSS::theSparseMatrix);
    long int& ldx(*_ldx);
    long int& ldy(*_ldy);
    int& blockSize(*_blockSize);

    for (int ib=0; ib<blockSize; ib++)
    {
        T* xvec = static_cast<T*>(x) + ldx*ib;
        T* yvec = static_cast<T*>(y) + ldy*ib;
        SSS::theSparseMatrix->DoMVMultiplication(primme->m,primme->n,xvec,yvec,*transpose);
    }
    *ierr = 0;
}

template <class T> void DenseMatvec(void *x, PRIMME_INT *_ldx, void *y, PRIMME_INT *_ldy, int *_blockSize,int* transpose,  primme_svds_params *primme, int *ierr)
{
    typedef SparseSVDSolver<T> SSS;
    assert(SSS::theDenseMatrix);
    long int& ldx(*_ldx);
    long int& ldy(*_ldy);
    int& blockSize(*_blockSize);
    int M=primme->m;
    int N=primme->n;
    if (*transpose==0)
        for (int ib=0; ib<blockSize; ib++)
        {
            T* xvec = static_cast<T*>(x) + ldx*ib;
            T* yvec = static_cast<T*>(y) + ldy*ib;
            for (int i=1;i<=M;i++)
            {
                yvec[i-1]=0.0;
                for (int j=1;j<=N;j++)
                    yvec[i-1]+=(*SSS::theDenseMatrix)(i,j)*xvec[j-1];
            }
        }
    else
        for (int ib=0; ib<blockSize; ib++)
        {
            T* xvec = static_cast<T*>(x) + ldx*ib;
            T* yvec = static_cast<T*>(y) + ldy*ib;
            for (int i=1;i<=N;i++)
            {
                yvec[i-1]=0.0;
                for (int j=1;j<=M;j++)
                    yvec[i-1]+=conj((*SSS::theDenseMatrix)(j,i))*xvec[j-1];
            }
        }
    *ierr = 0;

}

template <class T>
void ClientMatvec(void *x, PRIMME_INT *_ldx, void *y, PRIMME_INT *_ldy, int *_blockSize, int* transpose, primme_svds_params *primme, int *ierr)
{
    assert(SparseSVDSolverClient<T>::theClient); //|Psi>
    long int& ldx(*_ldx);
    long int& ldy(*_ldy);
    int& blockSize(*_blockSize);

    for (int ib=0; ib<blockSize; ib++)
    {
        const T* xvec = static_cast<const T*>(x) + ldx*ib;
              T* yvec = static_cast<      T*>(y) + ldy*ib;
        SparseSVDSolverClient<T>::theClient->DoMatVecContraction(primme->n,xvec,yvec,*transpose);
    }
    *ierr = 0;

}


//
//  Make template instances
//
template class PrimeSVDSolver<std::complex<double> >;
template class PrimeSVDSolver<double>;

