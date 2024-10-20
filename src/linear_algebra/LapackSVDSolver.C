#include "oml/numeric/LapackSVDSolver.H"
#include "oml/matrix.h"
#include "oml/diagonalmatrix.h"
#include "oml/vector.h"
#include <complex>
//
// See http://www.netlib.org/lapack/explore-html/d1/d7e/group__double_g_esing_ga84fdf22a62b12ff364621e4713ce02f2.html
// for detailed docs
// you also need to add -llapack to the link command
typedef std::complex<double> dcmplx;

extern"C" {
void dgesvd_(char* JOBU,char* JOBVT,int* M,int* N,double* A,int* LDA, double* S,double* U,
int* LDU,double* VT,int* LDVT,double* WORK,int* LWORK,int* INFO);
void zgesvd_(char* JOBU,char* JOBVT,int* M,int* N,dcmplx* A,int* LDA, double* S,dcmplx* U,
int* LDU,dcmplx* VT,int* LDVT,dcmplx* WORK,int* LWORK,double* RWORK, int* INFO);
}


template <class T> void gesvd(char* JOBU,char* JOBV,int* M,int* N,T* A,int* LDA, double* S,T* U,
int* LDU,T* VT,int* LDVT,T* WORK,int* LWORK,int* INFO);

template <> void gesvd<double> (char* JOBU,char* JOBV,int* M,int* N,double* A,int* LDA, double* S,double* U,
int* LDU,double* VT,int* LDVT,double* WORK,int* LWORK,int* INFO)
{
    dgesvd_(JOBU,JOBV,M,N,A,LDA,S,U,LDU,VT,LDVT,WORK,LWORK,INFO); //double
}
template <> void gesvd<dcmplx> (char* JOBU,char* JOBV,int* M,int* N,dcmplx* A,int* LDA, double* S,dcmplx* U,
int* LDU,dcmplx* VT,int* LDVT,dcmplx* WORK,int* LWORK,int* INFO)
{
    Vector<double> rwork(5*Min(*M,*N));
    zgesvd_(JOBU,JOBV,M,N,A,LDA,S,U,LDU,VT,LDVT,WORK,LWORK,&rwork(1),INFO); //complex<double>
}

namespace oml {
template <class T> typename LapackSVDSolver<T>::UsVType
LapackSVDSolver<T>::SolveAll(const MatrixT& A)
{
    int mn=Min(A.GetNumRows(),A.GetNumCols());
    return Solve(A,mn);
}

template <class T> typename LapackSVDSolver<T>::UsVType
LapackSVDSolver<T>::Solve(const MatrixT& A,size_t NumSingularValues)
{
    assert(NumSingularValues>0);
    assert(NumSingularValues<=Min(A.GetNumRows(),A.GetNumCols()));
    assert(A.GetLimits().Row.Low==1);
    assert(A.GetLimits().Col.Low==1);
    assert(!isnan(A));
    int M=A.GetNumRows(),N=A.GetNumCols(),mn=Min(M,N);

    //
    //  Dicey deciding how much work space lapack needs. more is faster
    //
    int info=0,lwork=-1;
    Vector<double> s(mn);
    Vector<T> work(1);
    Matrix<T> U(A),VT(N,N);
    char jobu='O',jobv='A';
    //
    //  Initial call to see how much work space is needed
    //
    gesvd<T>(&jobu,&jobv,&M,&N,&U(1,1),&M,&s(1),0,&M,&VT(1,1),&N,&work(1),&lwork,&info);
    lwork=real(work(1));
    work.SetLimits(lwork);
    //
    //  Now do the actual SVD work
    //
    gesvd<T>(&jobu,&jobv,&M,&N,&U(1,1),&M,&s(1),0,&M,&VT(1,1),&N,&work(1),&lwork,&info);
    if (info!=0)
    {
        std::cerr << "Warning: LapackSVDSolver info=" << info << std::endl;
//        std::cout << std::fixed << "A=" << A << std::endl;
//        std::cout << "s=" << s << std::endl;
    }
    assert(info==0);
    //
    //  Now fix up the matrix limits
    //
    if (NumSingularValues < static_cast<unsigned>(mn))
    {
        s.SetLimits(NumSingularValues,true);
        mn=NumSingularValues;
    }
    U .SetLimits(M,mn,true); //Throw away last N-mn columns
    VT.SetLimits(mn,N,true); // Throw away last N-mn rows
    DiagonalMatrix<double> ds(s);
    return std::make_tuple(std::move(U),std::move(ds),std::move(VT));
}

template class LapackSVDSolver<double>;
template class LapackSVDSolver<dcmplx>;

} //namespace oml
