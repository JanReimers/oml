#include "oml/numeric/LapackEigenSolver.H"
#include "oml/matrix.h"
#include "oml/vector.h"
#include "oml/fakedouble.h"
#include <complex>
#include <iostream>
//
// See http://www.netlib.org/lapack/explore-html/d1/d7e/group__double_g_esing_ga84fdf22a62b12ff364621e4713ce02f2.html
// for detailed docs
// handy search tool: https://software.intel.com/content/www/us/en/develop/articles/intel-mkl-function-finding-advisor.html
// you also need to add -llapack to the link command
typedef std::complex<double> dcmplx;

extern"C" {

void dsyevx_(char* JOBZ,char* RANGE,char* UPLO,int* N,double* A,int* LDA,double* VL,double* VU,int* IL,int* IU,double* ABSTOL,int* M,double* W,double* Z,int* LDZ,double* WORK,int* LWORK,              int* IWORK,int* IFAIL,int* INFO);
void zheevx_(char* JOBZ,char* RANGE,char* UPLO,int* N,dcmplx* A,int* LDA,double* VL,double* VU,int* IL,int* IU,double* ABSTOL,int* M,double* W,dcmplx* Z,int* LDZ,dcmplx* WORK,int* LWORK,double* RWORK,int* IWORK,int* IFAIL,int* INFO);
void dgeev_ (char* JOBVL,char* JOBVR,int* N,double* A,int* LDA,double* WR,double* WI,double* VL,int* LDVL,double* VR,int* LDVR,double* WORK,int* LWORK,              int* INFO);
void zgeev_ (char* JOBVL,char* JOBVR,int* N,dcmplx* A,int* LDA,dcmplx* W            ,dcmplx* VL,int* LDVL,dcmplx* VR,int* LDVR,dcmplx* WORK,int* LWORK,double* RWORK,int* INFO);
}
// Symmetric
template <class T> void evx  (char* JOBZ,char* RANGE,char* UPLO,int* N,T     * A,int* LDA,double* VL,double* VU,int* IL,int* IU,double* ABSTOL,int* M,double* W,T     * Z,int* LDZ,T     * WORK,int* LWORK,int* IWORK,int* IFAIL,int* INFO);
template <> void evx<double> (char* JOBZ,char* RANGE,char* UPLO,int* N,double* A,int* LDA,double* VL,double* VU,int* IL,int* IU,double* ABSTOL,int* M,double* W,double* Z,int* LDZ,double* WORK,int* LWORK,int* IWORK,int* IFAIL,int* INFO)
{
    dsyevx_(JOBZ,RANGE,UPLO,N,A,LDA,VL,VU,IL,IU,ABSTOL,M,W,Z,LDZ,WORK,LWORK,IWORK,IFAIL,INFO); //double
}
template <> void evx<dcmplx> (char* JOBZ,char* RANGE,char* UPLO,int* N,dcmplx* A,int* LDA,double* VL,double* VU,int* IL,int* IU,double* ABSTOL,int* M,double* W,dcmplx* Z,int* LDZ,dcmplx* WORK,int* LWORK,int* IWORK,int* IFAIL,int* INFO)
{
    Vector<double> rwork(7*(*N));
    zheevx_(JOBZ,RANGE,UPLO,N,A,LDA,VL,VU,IL,IU,ABSTOL,M,W,Z,LDZ,WORK,LWORK,&rwork(1),IWORK,IFAIL,INFO); //complex<double>
}
// Non symmetric
template <class T> void ev (char* JOBVL,char* JOBVR,int* N,T     * A,int* LDA,T     * WR,T     * WI,T     * VL,int* LDVL,T     * VR,int* LDVR,T     * WORK,int* LWORK,int* INFO);
template <> void ev<double>(char* JOBVL,char* JOBVR,int* N,double* A,int* LDA,double* WR,double* WI,double* VL,int* LDVL,double* VR,int* LDVR,double* WORK,int* LWORK,int* INFO)
{
    dgeev_ (JOBVL,JOBVR,N,A,LDA,WR,WI,VL,LDVL,VR,LDVR,WORK,LWORK,      INFO);
}
template <> void ev<dcmplx>(char* JOBVL,char* JOBVR,int* N,dcmplx* A,int* LDA,dcmplx* W,dcmplx*    ,dcmplx* VL,int* LDVL,dcmplx* VR,int* LDVR,dcmplx* WORK,int* LWORK,int* INFO)
{
    Vector<double> rwork(7*(*N));
    zgeev_ (JOBVL,JOBVR,N,A,LDA,W    ,VL,LDVL,VR,LDVR,WORK,LWORK,&rwork(1),INFO);
}

namespace oml {

template <class T> typename LapackEigenSolver<T>::UdType
LapackEigenSolver<T>::SolveAll(const MatrixT& A,double eps)
{
    int N=A.GetNumRows();
    if (IsDiagonal(A,N*eps))
    {
        Matrix<T> U(N,N);
        Unit(U);
        Vector<T> e=A.GetDiagonal();
        return std::make_tuple(U,real(e));
    }
    else
        return Solve(A,eps,N);
}



template <class T> typename LapackEigenSolver<T>::UdType
LapackEigenSolver<T>::Solve(const MatrixT& A,double eps, int NumEigenValues)
{
    assert(A.GetLimits().Row.Low==1);
    assert(A.GetLimits().Col.Low==1);
    int N=A.GetNumRows();
    assert(static_cast<unsigned>(N)==A.GetNumCols());
    assert(NumEigenValues<=N);
    double epsH=eps == 0.0 ? 1e-13 : eps;
    if (!IsHermitian(A,epsH))
    {
        std::cerr << "Non-hermitian matrix, eps=" << eps << std::endl;
        double delta=Max(fabs(A-Transpose(conj(A))));
        std::cerr << "Max(fabs(A-Adagger))=" << delta << std::endl;
    }
    assert(IsHermitian(A,epsH));
    //
    //  Dicey deciding how much work space lapack needs. more is faster
    //
    int info=0,lwork=-1,IL=1,IU=NumEigenValues;
    double VL,VU;
    Vector<double> W(N);
    Vector<T> work(1);
    Vector<int> iwork(5*N),ifail(N);
    Matrix<T> U(N,NumEigenValues),Alower(A);
    char jobz='V',range='I',uplo='L';
    //
    //  Initial call to see how much work space is needed
    //
    evx<T>(&jobz,&range,&uplo,&N,&Alower(1,1),&N,&VL,&VU,&IL,&IU,&eps,&NumEigenValues,&W(1),&U(1,1),&N, &work(1),&lwork,&iwork(1),&ifail(1),&info);
    lwork=real(work(1));
    work.SetLimits(lwork);
    //
    //  Now do the actual SVD work
    //
    evx<T>(&jobz,&range,&uplo,&N,&Alower(1,1),&N,&VL,&VU,&IL,&IU,&eps,&NumEigenValues,&W(1),&U(1,1),&N, &work(1),&lwork,&iwork(1),&ifail(1),&info);
    if (info!=0)
    {
        std::cerr << "Warning: LapackEigenSolver info=" << info << std::endl;
//        std::cout << std::fixed << "A=" << A << std::endl;
//        std::cout << "W=" << W << std::endl;
    }
//    assert(info==0);
    //
    //  Now fix up the matrix limits
    W.SetLimits(NumEigenValues,true);
    return std::make_tuple(std::move(U),std::move(W));
}

template <class T> typename LapackEigenSolver<T>::UdTypeN
LapackEigenSolver<T>::SolveAllRightNonSym(const MatrixT& A,double eps)
{
    size_t N=A.GetNumRows();
    return SolveRightNonSym(A,eps,N);
}

template <> typename LapackEigenSolver<double>::UdTypeN
LapackEigenSolver<double>::SolveRightNonSym(const MatrixT& A,double eps, int NumEigenValues)
{
    assert(A.GetLimits().Row.Low==1);
    assert(A.GetLimits().Col.Low==1);
    int N=A.GetNumRows();
    assert(static_cast<unsigned>(N)==A.GetNumCols());
    if (NumEigenValues<N)
        std::cerr << "Warning: Lapack does not support subset of eigen values for non-symmtric matrcies" << std::endl;
    //
    //  Dicey deciding how much work space lapack needs. more is faster
    //
    int info=0,lwork=-1;
    Vector<double> WR(N),WI(N);
    Vector<double> work(1);
    Vector<int> iwork(5*N),ifail(N);
    Matrix<double> VL(N,N),VR(N,N),Acopy(A);
    char jobvl='N',jobvr='V';
    //
    //  Initial call to see how much work space is needed
    //
    dgeev_(&jobvl,&jobvr,&N,&Acopy(1,1),&N,&WR(1),&WI(1),&VL(1,1),&N,&VR(1,1),&N, &work(1),&lwork,&info);
    lwork=real(work(1));
    work.SetLimits(lwork);
    //
    //  Now do the actual Eigen decomposition work
    //
    dgeev_(&jobvl,&jobvr,&N,&Acopy(1,1),&N,&WR(1),&WI(1),&VL(1,1),&N,&VR(1,1),&N, &work(1),&lwork,&info);
    if (info!=0)
        std::cerr << "Warning: LapackEigenSolver info=" << info << std::endl;
    assert(info==0);
    //
    //  Unpack the eigen pairs
    //
    Vector<dcmplx> W(N);
    Matrix<dcmplx> V(N,N);
    for (int j=1;j<=N;j++)
    {
        if (fabs(WI(j))<eps)
        {
            W(j)=dcmplx(WR(j),WI(j));
            for (int i=1;i<=N;i++)
                V(i,j)=dcmplx(VR(i,j),0.0);
        }
        else
        {
            W(j  )=dcmplx(WR(j  ),WI(j  ));
            W(j+1)=dcmplx(WR(j+1),WI(j+1));
            for (int i=1;i<=N;i++)
            {
                V(i,j  )=dcmplx(VR(i,j), VR(i,j+1));
                V(i,j+1)=dcmplx(VR(i,j),-VR(i,j+1));
            }
            j++;
        }
    }
    return std::make_tuple(std::move(V),std::move(W));
}

template <> typename LapackEigenSolver<dcmplx>::UdTypeN
LapackEigenSolver<dcmplx>::SolveRightNonSym(const MatrixT& A,double eps, int NumEigenValues)
{
    assert(A.GetLimits().Row.Low==1);
    assert(A.GetLimits().Col.Low==1);
    int N=A.GetNumRows();
    assert(static_cast<unsigned>(N)==A.GetNumCols());
    if (NumEigenValues<N)
        std::cerr << "Warning: Lapack does not support subset of eigen values for non-symmtric matrcies" << std::endl;
    //
    //  Dicey deciding how much work space lapack needs. more is faster
    //
    int info=0,lwork=-1;
    Vector<dcmplx> W(N),work(1);
    Vector<double> rwork(2*N);
    Matrix<dcmplx> VL(N,N),VR(N,N),Acopy(A);
    char jobvl='N',jobvr='V';
    //
    //  Initial call to see how much work space is needed
    //
    zgeev_(&jobvl,&jobvr,&N,&Acopy(1,1),&N,&W(1),&VL(1,1),&N,&VR(1,1),&N, &work(1),&lwork,&rwork(1),&info);
    lwork=real(work(1));
    work.SetLimits(lwork);
    //
    //  Now do the actual SVD work
    //
    zgeev_(&jobvl,&jobvr,&N,&Acopy(1,1),&N,&W(1),&VL(1,1),&N,&VR(1,1),&N, &work(1),&lwork,&rwork(1),&info);
    assert(info==0);
    return std::make_tuple(std::move(VR),std::move(W));
}

//
//  Make template instances
//
template class LapackEigenSolver<std::complex<double> >;
template class LapackEigenSolver<double>;

} //namespace oml
