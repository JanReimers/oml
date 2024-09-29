
#include "ArpackEigenSolver.H"
#include "Containers/SparseMatrix.H"

#include "arpack/arpackdef.h"
#include "oml/matrix.h"
#include <tuple>
#include <iostream>

using std::cout;
using std::endl;

typedef std::complex<double> dcmplx;

//
//  Hand modified version of arpack.h with the horrible C-style __complex__ crap removed.
//  This code uses the fork arpack-ng https://github.com/opencollab/arpack-ng
//  These routines are very hard ot use becuase they force the user to implement the interation loop
//  They dont seem to undestand how to pass function names for matvec operations like lapack and primme do.
//  See here for docs: https://www.caam.rice.edu/software/ARPACK/UG/node138.html
//                     https://www.caam.rice.edu/software/ARPACK/UG/node44.html
//                      https://scm.cs.kuleuven.be/scm/svn/numerics_software/ARPACK/SRC/dneupd.f//
//
//  INFO!= error codes:
//c  INFO    Integer.  (INPUT/OUTPUT)
//c          If INFO .EQ. 0, a randomly initial residual vector is used.
//c          If INFO .NE. 0, RESID contains the initial residual vector,
//c                          possibly from a previous run.
//c          Error flag on output.
//c          =  0: Normal exit.
//c          =  1: Maximum number of iterations taken.
//c                All possible eigenvalues of OP has been found. IPARAM(5)
//c                returns the number of wanted converged Ritz values.
//c          =  2: No longer an informational error. Deprecated starting
//c                with release 2 of ARPACK.
//c          =  3: No shifts could be applied during a cycle of the
//c                Implicitly restarted Arnoldi iteration. One possibility
//c                is to increase the size of NCV relative to NEV.
//c                See remark 4 below.
//c          = -1: N must be positive.
//c          = -2: NEV must be positive.
//c          = -3: NCV-NEV >= 2 and less than or equal to N.
//c          = -4: The maximum number of Arnoldi update iteration
//c                must be greater than zero.
//c          = -5: WHICH must be one of 'LM', 'SM', 'LR', 'SR', 'LI', 'SI'
//c          = -6: BMAT must be one of 'I' or 'G'.
//c          = -7: Length of private work array is not sufficient.
//c          = -8: Error return from LAPACK eigenvalue calculation;
//c          = -9: Starting vector is zero.
//c          = -10: IPARAM(7) must be 1,2,3,4.
//c          = -11: IPARAM(7) = 1 and BMAT = 'G' are incompatable.
//c          = -12: IPARAM(1) must be equal to 0 or 1.
//c          = -9999: Could not build an Arnoldi factorization.
//c                   IPARAM(5) returns the size of the current Arnoldi
//c                   factorization.
//c
extern "C"
{
    // Non symmetric complex single
    void cnaupd_c(a_int* ido, char const* bmat, a_int n, char const* which, a_int nev, float  tol, std::complex<float>* resid, a_int ncv, std::complex<float>* v, a_int ldv, a_int* iparam, a_int* ipntr, std::complex<float>* workd, std::complex<float>* workl, a_int lworkl, float* rwork, a_int* info);
    // Non symmetric double
    void dnaupd_c(a_int* ido, char const* bmat, a_int n, char const* which, a_int nev, double tol, double* resid, a_int ncv, double* v, a_int ldv, a_int* iparam, a_int* ipntr, double* workd, double* workl, a_int lworkl, a_int* info);
    // Symmetric double
    void dsaupd_c(a_int* ido, char const* bmat, a_int n, char const* which, a_int nev, double tol, double* resid, a_int ncv, double* v, a_int ldv, a_int* iparam, a_int* ipntr, double* workd, double* workl, a_int lworkl, a_int* info);
    // Non symmetric single
    void snaupd_c(a_int* ido, char const* bmat, a_int n, char const* which, a_int nev, float  tol, float* resid, a_int ncv, float* v, a_int ldv, a_int* iparam, a_int* ipntr, float* workd, float* workl, a_int lworkl, a_int* info);
    // Symmetric single
    void ssaupd_c(a_int* ido, char const* bmat, a_int n, char const* which, a_int nev, float  tol, float* resid, a_int ncv, float* v, a_int ldv, a_int* iparam, a_int* ipntr, float* workd, float* workl, a_int lworkl, a_int* info);
    // Non symmetric complex double
    void znaupd_c(a_int* ido, char const* bmat, a_int n, char const* which, a_int nev, double tol, dcmplx* resid, a_int ncv, dcmplx* v, a_int ldv, a_int* iparam, a_int* ipntr, dcmplx* workd, dcmplx* workl, a_int lworkl, double* rwork, a_int* info);

    // Call these to extract eigen vectors after the iterations are complete
    void cneupd_c(bool rvec, char const* howmny, a_int const* select, std::complex<float>*  d, std::complex<float>* z, a_int ldz, std::complex<float> sigma, std::complex<float>* workev, char const* bmat, a_int n, char const* which, a_int nev, float tol, std::complex<float>* resid, a_int ncv, std::complex<float>* v, a_int ldv, a_int* iparam, a_int* ipntr, std::complex<float>* workd, std::complex<float>* workl, a_int lworkl, float* rwork, a_int* info);
    void dneupd_c(bool rvec, char const* howmny, a_int const* select, double* dr, double* di, double* z, a_int ldz, double sigmar, double sigmai, double * workev, char const* bmat, a_int n, char const* which, a_int nev, double tol, double* resid, a_int ncv, double* v, a_int ldv, a_int* iparam, a_int* ipntr, double* workd, double* workl, a_int lworkl, a_int* info);
    void dseupd_c(bool rvec, char const* howmny, a_int const* select, double* d , double*  z, a_int ldz, double sigma, char const* bmat, a_int n, char const* which, a_int nev, double tol, double* resid, a_int ncv, double* v, a_int ldv, a_int* iparam, a_int* ipntr, double* workd, double* workl, a_int lworkl, a_int* info);
    void sneupd_c(bool rvec, char const* howmny, a_int const* select, float*  dr, float*  di, float* z, a_int ldz, float sigmar, float sigmai, float * workev, char const* bmat, a_int n, char const* which, a_int nev, float tol, float* resid, a_int ncv, float* v, a_int ldv, a_int* iparam, a_int* ipntr, float* workd, float* workl, a_int lworkl, a_int* info);
    void sseupd_c(bool rvec, char const* howmny, a_int const* select, float*  d , float*   z, a_int ldz, float sigma, char const* bmat, a_int n, char const* which, a_int nev, float tol, float* resid, a_int ncv, float* v, a_int ldv, a_int* iparam, a_int* ipntr, float* workd, float* workl, a_int lworkl, a_int* info);
    void zneupd_c(bool rvec, char const* howmny, a_int const* select, dcmplx* d , dcmplx*  z, a_int ldz, dcmplx sigma, dcmplx* workev, char const* bmat, a_int n, char const* which, a_int nev, double tol, dcmplx* resid, a_int ncv, dcmplx* v, a_int ldv, a_int* iparam, a_int* ipntr, dcmplx* workd, dcmplx* workl, a_int lworkl, double* rwork, a_int* info);
}

//----------------------------------------------------------------------------------------
//
//  Private solver routines routines
//
template <> typename ArpackEigenSolver<dcmplx>::UdTypeN
ArpackEigenSolver<dcmplx>::SolveG(MatvecT matvec,int N, int Nev,double eps)
{
    assert(Nev>0);
    assert(Nev<=N-2);

    int IDO=0,INFO=0;
    int Ncv=5*Nev; //  *** THis has a huge effect on convergence, bigger is better.
    if (Ncv>N) Ncv=N;
    int Lworkl=3*Ncv*Ncv + 5*Ncv;

    int MaxIter=1000;

    Vector<dcmplx>  residuals(N),Workd(3*N),Workl(Lworkl);
    Matrix<dcmplx>  V(N,Ncv);
    Vector<double> rwork(Ncv);
    Vector<int> iParam(11);
    iParam(1)=1; //ISHIFT
    iParam(3)=MaxIter;
    iParam(4)=1; //NB only 1 works
    iParam(7)=1; //Mode
    int iPntr[14];
    char arI='I';
    char which[] = "LR";

    // Arnaldi iteration loop
    do
    {
        znaupd_c(&IDO,&arI,N,which,Nev,eps,&residuals(1),Ncv,&V(1,1)
        ,N,&iParam(1),iPntr,&Workd(1),&Workl(1),Lworkl,&rwork(1),&INFO);
//        cout << "IDO=" << IDO << endl;
        if (IDO==-1 || IDO==1)
            matvec(N,&Workd(iPntr[0]),&Workd(iPntr[1]));
        else
            break;

    } while(true);
    if (INFO!=0)
    {
        cout << "Info=" << INFO << endl;
    }
    //cout << "Info=" << INFO << endl;
    //cout << "nIter=" << iParam(3) << endl;
    assert(INFO==0);
    //
    //  Post processing to the eigen vectors.
    //
    Vector<int> select(Ncv);
    Vector<dcmplx> D(Nev+1),Workev(2*Ncv);
    Matrix<dcmplx> U(N,Nev);
    dcmplx sigma(0);
    zneupd_c(true, "All", &select(1), &D(1), &U(1,1),N,
        sigma, &Workev(1), &arI,N,which, Nev, eps,
          &residuals(1), Ncv, &V(1,1),N, &iParam(1), iPntr, &Workd(1), &Workl(1),
          Lworkl, &rwork(1), &INFO );

    theDenseMatrix=nullptr; //Make sure these don't accidentally get reused
    theSparseMatrix=nullptr;
    ClientT::theClient=nullptr;
//    cout << "D(Nev+1)=" << D(Nev+1) << endl;
    D.SetLimits(Nev,true);
    return make_tuple(std::move(U),std::move(D));
}



template <> typename ArpackEigenSolver<double>::UdTypeN
ArpackEigenSolver<double>::SolveG(MatvecT matvec,int N, int Nev,double eps)
{
    assert(Nev>0);
    assert(Nev<=N-2);

    int IDO=0,INFO=0;
    int Ncv=2*Nev+1; //*** THis has a huge effect on convergence, bigger is better.
    if (Ncv>=N) Ncv=N;
    //cout << "Nev Ncv N = " << Nev << " " << Ncv << " " << N << endl;
    int Lworkl=3*Ncv*Ncv + 6*Ncv;

    int MaxIter=1000;

    Vector<double>  residuals(N);
    Vector<double>  Workd(3*N);
    Vector<double>  Workl(Lworkl);
    Matrix<double>  V(N,Ncv);
    Vector<double> rwork(Ncv);
    Vector<int> iParam(11);
    iParam(1)=1; //ISHIFT
    iParam(3)=MaxIter;
    iParam(4)=1; //NB only 1 works
    iParam(7)=1; //Mode
    int iPntr[14];
    char arI='I';
    char which[] = "LR";

    // Arnaldi iteration loop
    do
    {
        dnaupd_c(&IDO,&arI,N,which,Nev,eps,&residuals(1),Ncv,&V(1,1)
        ,N,&iParam(1),iPntr,&Workd(1),&Workl(1),Lworkl,&INFO);
//        cout << "IDO=" << IDO << endl;
        if (IDO==-1 || IDO==1)
            matvec(N,&Workd(iPntr[0]),&Workd(iPntr[1]));
        else
            break;

    } while(true);
    if (INFO!=0)
    {
        cout << "Info=" << INFO << endl;
    }
//    cout << "nIter=" << iParam(3) << endl;
    assert(INFO==0);
    //
    //  Post processing to the eigen vectors.
    //
    Vector<int> select(Ncv);
    Vector<double> DR(Nev+1),DI(Nev+1),Workev(3*Ncv);
    Fill(DI,0.0);
    double sigmar(0),sigmai(0);
    char how_many('A');
    dneupd_c(true, &how_many, &select(1), &DR(1),&DI(1), &V(1,1),N,
        sigmar,sigmai, &Workev(1), &arI,N,which, Nev, eps,
          &residuals(1), Ncv, &V(1,1),N, &iParam(1), iPntr, &Workd(1), &Workl(1),
          Lworkl, &INFO );

    if (fabs(DI(Nev+1))>eps)
    {
        std::cerr << "warning Incrementing Nev to avoid conjugate eigen value pair cut off" << std::endl;
        Nev++; //avoid cutting off a conj pair.
    }
    Vector<dcmplx> D(Nev);
    Matrix<dcmplx> UC(N,Nev);
    for (int j=1;j<=Nev;j++)
    {
        if (fabs(DI(j))>eps && j==Nev)
            std::cerr << "warning conjugate eigen value pair is cut off, increase Nev by one" << std::endl;
        if (fabs(DI(j))<eps || j==Nev)
        {
            D(j)=dcmplx(DR(j),0.0);
            for (int i=1;i<=N;i++)
                UC(i,j  )=dcmplx(V(i,j),0.0);
        }
        else
        {
            D(j  )=dcmplx(DR(j  ),DI(j  ));
            D(j+1)=dcmplx(DR(j+1),DI(j+1));
            for (int i=1;i<=N;i++)
            {
                UC(i,j  )=dcmplx(V(i,j), V(i,j+1));
                UC(i,j+1)=dcmplx(V(i,j),-V(i,j+1));
            }
            j++;
        }
    }
    theDenseMatrix=nullptr; //Make sure these don't accidentally get reused
    theSparseMatrix=nullptr;
    ClientT::theClient=nullptr;
    return make_tuple(std::move(UC),std::move(D));
}

template <> typename ArpackEigenSolver<double>::UdType
ArpackEigenSolver<double>::SolveSym(MatvecT matvec,int N, int Nev,double eps)
{
    assert(Nev>0);
    assert(Nev<=N-2);

    int IDO=0,INFO=0;
    int Ncv=2*Nev+1; //*** THis has a huge effect on convergence, bigger is better.
    if (Ncv>=N) Ncv=N;
    //cout << "Nev Ncv N = " << Nev << " " << Ncv << " " << N << endl;
    int Lworkl=Ncv*Ncv + 8*Ncv;

    int MaxIter=1000;

    Vector<double>  residuals(N);
    Vector<double>  Workd(3*N);
    Vector<double>  Workl(Lworkl);
    Matrix<double>  V(N,Ncv);
    Vector<double> rwork(Ncv);
    Vector<int> iParam(11);
    iParam(1)=1; //ISHIFT
    iParam(3)=MaxIter;
    iParam(4)=1; //NB only 1 works
    iParam(7)=1; //Mode
    int iPntr[14];
    char arI='I';

    // Arnaldi iteration loop
    do
    {
      dsaupd_c(&IDO,&arI,N,"LA",Nev,eps,&residuals(1),Ncv,&V(1,1),N,&iParam(1),iPntr,&Workd(1),&Workl(1),Lworkl,&INFO);
//        cout << "IDO=" << IDO << endl;
        if (IDO==-1 || IDO==1)
            matvec(N,&Workd(iPntr[0]),&Workd(iPntr[1]));
        else
            break;

    } while(true);
    if (INFO!=0)
    {
        cout << "Info=" << INFO << endl;
        cout << "nIter=" << iParam(3) << endl;
    }
    assert(INFO==0);
    //
    //  Post processing to the eigen vectors.
    //
    Vector<int> select(Ncv);
    Vector<double> D(Nev+1),Workev(3*Ncv);
    double sigma(0);
    char how_many('A');
    dseupd_c(true, &how_many, &select(1), &D(1), &V(1,1),N, sigma, &arI,N,"LA", Nev, eps, &residuals(1), Ncv, &V(1,1),N, &iParam(1), iPntr, &Workd(1), &Workl(1),          Lworkl, &INFO );

    D.SetLimits(Nev,true);
    V.SetLimits(N,Nev,true);
    theDenseMatrix=nullptr; //Make sure these don't accidentally get reused
    theSparseMatrix=nullptr;
    ClientT::theClient=nullptr;
    return std::make_tuple(std::move(V),std::move(D));
}



template <class T> void matvecDense(int N, const T * x, T * y)
{
    typedef ArpackEigenSolver<T> ArpackT;
    assert(ArpackT::theDenseMatrix);
    assert(ArpackT::theDenseMatrix->GetNumRows()==N);
    assert(ArpackT::theDenseMatrix->GetNumCols()==N);
    for (int i=1;i<=N;i++)
    {
        y[i-1]=T (0.0);
        for (int j=1;j<=N;j++)
            y[i-1]+=(*ArpackT::theDenseMatrix)(i,j)*x[j-1];
    }
}

template <class T> void matvecSparse(int N, const T * x, T * y)
{
    typedef ArpackEigenSolver<T> ArpackT;
    assert(ArpackT::theSparseMatrix);
    assert(ArpackT::theSparseMatrix->GetNumRows()==N);
    assert(ArpackT::theSparseMatrix->GetNumCols()==N);
    ArpackT::theSparseMatrix->DoMVMultiplication(N,x,y);
}

template <class T> void matvecClient(int N, const T * x, T * y)
{
    typedef typename ArpackEigenSolver<T>::ClientT ClientT;
    assert(ClientT::theClient);
    assert(ClientT::theClient->GetSize()==N);
    ClientT::theClient->DoMatVecContraction(N,x,y);
}

//----------------------------------------------------------------------------------------
//
//  Public interface routines
//
template <class T> typename ArpackEigenSolver<T>::UdType
ArpackEigenSolver<T>::Solve(const MatrixT& A,double eps, int Nev)
{
    // For complex Hermitian we are supposed to use the non-sym solver
    assert(A.GetLimits().Row.Low==1);
    assert(A.GetLimits().Col.Low==1);
    assert(A.GetNumRows()==A.GetNumCols());
    assert(Nev<A.GetNumRows());
    theDenseMatrix=&A;
    auto [Uc,dc]=SolveG(matvecDense,A.GetNumRows(),Nev,eps);
    Vector<double> d=real(dc);
    return std::make_tuple(std::move(Uc),std::move(d));
}

template <> typename ArpackEigenSolver<double>::UdType
ArpackEigenSolver<double>::Solve(const MatrixT& A,double eps, int Nev)
{
    assert(A.GetLimits().Row.Low==1);
    assert(A.GetLimits().Col.Low==1);
    assert(A.GetNumRows()==A.GetNumCols());
    assert(Nev<A.GetNumRows());
    assert(IsSymmetric(A,eps));
    theDenseMatrix=&A;
    return SolveSym(matvecDense,A.GetNumRows(),Nev,eps);
}

template <class T> typename ArpackEigenSolver<T>::UdType  ArpackEigenSolver<T>::SolveAll      (const MatrixT& A,double eps)
{
    int N=A.GetNumRows();
    return Solve(A,eps,N-2);
}



template <class T> typename ArpackEigenSolver<T>::UdType
ArpackEigenSolver<T>::Solve(const SparseMatrixT& A,double eps, int Nev)
{
    // For complex Hermitian we are supposed to use the non-sym solver
    assert(A.GetLimits().Row.Low==1);
    assert(A.GetLimits().Col.Low==1);
    assert(A.GetNumRows()==A.GetNumCols());
    assert(Nev<A.GetNumRows());
    theSparseMatrix=&A;
    auto [Uc,dc]=SolveG(matvecSparse,A.GetNumRows(),Nev,eps);
    Vector<double> d=real(dc);
    return std::make_tuple(std::move(Uc),std::move(d));
}

template <> typename ArpackEigenSolver<double>::UdType
ArpackEigenSolver<double>::Solve(const SparseMatrixT& A,double eps, int Nev)
{
    assert(A.GetLimits().Row.Low==1);
    assert(A.GetLimits().Col.Low==1);
    assert(A.GetNumRows()==A.GetNumCols());
    assert(Nev<A.GetNumRows());
    theSparseMatrix=&A;
    return SolveSym(matvecSparse,A.GetNumRows(),Nev,eps);
}


template <class T> typename ArpackEigenSolver<T>::UdTypeN
ArpackEigenSolver<T>::SolveRightNonSym(const MatrixT& A,double eps, int Nev)
{
    assert(A.GetLimits().Row.Low==1);
    assert(A.GetLimits().Col.Low==1);
    assert(A.GetNumRows()==A.GetNumCols());
    assert(Nev<A.GetNumRows());
    theDenseMatrix=&A;
    return SolveG(matvecDense,A.GetNumRows(),Nev,eps);
}


template <class T> typename ArpackEigenSolver<T>::UdTypeN ArpackEigenSolver<T>::SolveAllRightNonSym(const MatrixT& A,double eps)
{
    std::cerr << "ArpackEigenSolver does not support all Evs, doing N-2" << std::endl;
    int N=A.GetNumRows();
    return SolveRightNonSym(A,eps,N-2);
}

template <class T> typename ArpackEigenSolver<T>::UdTypeN
ArpackEigenSolver<T>::SolveRightNonSym(const SparseMatrixT& A,double eps, int Nev)
{
    assert(A.GetLimits().Row.Low==1);
    assert(A.GetLimits().Col.Low==1);
    assert(A.GetNumRows()==A.GetNumCols());
    assert(Nev<A.GetNumRows());
    theSparseMatrix=&A;
    return SolveG(matvecSparse,A.GetNumRows(),Nev,eps);
}

template <class T> typename ArpackEigenSolver<T>::UdType
ArpackEigenSolver<T>::Solve(const ClientT* client,double eps, int Nev)
{
    // For complex Hermitian we are supposed to use the non-sym solver
    ClientT::theClient=client;
    auto [Uc,dc]=SolveG(matvecClient,client->GetSize(),Nev,eps);
    Vector<double> d=real(dc);
    return std::make_tuple(std::move(Uc),std::move(d));
}

template <> typename ArpackEigenSolver<double>::UdType
ArpackEigenSolver<double>::Solve(const ClientT* client,double eps, int Nev)
{
    ClientT::theClient=client;
    return SolveSym(matvecClient,client->GetSize(),Nev,eps);
}


template <class T> typename ArpackEigenSolver<T>::UdTypeN
ArpackEigenSolver<T>::SolveRightNonSym(const ClientT* client,double eps, int Nev)
{
    ClientT::theClient=client;
    return SolveG(matvecClient,client->GetSize(),Nev,eps);
}


//
//  static variable. Kludge for getting the matrix into the MatVec routines.
//
template <class T> const SparseMatrix<T>* ArpackEigenSolver<T>::theSparseMatrix = 0;
template <class T> const       Matrix<T>* ArpackEigenSolver<T>::theDenseMatrix = 0;

template class ArpackEigenSolver<double>;
template class ArpackEigenSolver<dcmplx>;

