module;
#include <cassert>
#include <cmath>
#include <iostream>
#include <complex>

export module oml.Lapack:Linear;
import oml.Vector;
import oml.Matrix;
import oml.Solvers;

export namespace oml
{
//
//  Abstract interface for QR solvers
//
template <class T> class LapackLinearSolver
    : public virtual LinearSolver<T>
{
protected:
    using MatrixT=typename LinearSolver<T>::MatrixT;
    using VectorT=typename LinearSolver<T>::VectorT;
public:
    virtual ~LapackLinearSolver() {};
//
//  Non packed triangular systems, solve
//
    virtual VectorT SolveUpperTri(const MatrixT& A,const VectorT& b); // A is upper triangular
    virtual VectorT SolveLowerTri(const MatrixT& A,const VectorT& b); // A is lower triangular
    virtual VectorT Solve(const MatrixT& A,const VectorT& b); // A*x=b
//    virtual VectorT Solve(const VectorT& b,const MatrixT& A) {return LinearSolver<T>::Solve(b,A);}   // b=x*A

private:
    VectorT SolveTri(const MatrixT& A,const VectorT& b, char UL);
};

}

//
// See http://www.netlib.org/lapack/explore-html/da/dba/group__double_o_t_h_e_rcomputational_ga7068947990361e55177155d044435a5c.html
// for detailed docs
// you also need to add -llapack to the link command
typedef std::complex<double> dcmplx;

extern"C" {
void dtrtrs_(char* UPLO,char* TRANS, char* DIAG, int* N, int* NRHS, const double* A,int* LDA,double* B,int* LDB,int* INFO);
void ztrtrs_(char* UPLO,char* TRANS, char* DIAG, int* N, int* NRHS, const dcmplx* A,int* LDA,dcmplx* B,int* LDB,int* INFO);
void  dgesv_(int* N, int* NRHS,double* A, int* LDA, int* IPIV, const double* B, int* LDB, int* INFO);
void  zgesv_(int* N, int* NRHS,dcmplx* A, int* LDA, int* IPIV, const dcmplx* B, int* LDB, int* INFO);
}

template <class T> void xgesv(int* N, int* NRHS,T* A, int* LDA, int* IPIV, T* B, int* LDB, int* INFO);
template <> void xgesv<double>(int* N, int* NRHS,double* A, int* LDA, int* IPIV, double* B, int* LDB, int* INFO)
{
    dgesv_(N,NRHS,A,LDA,IPIV,B,LDB,INFO);
}
template <> void xgesv<dcmplx>(int* N, int* NRHS,dcmplx* A, int* LDA, int* IPIV, dcmplx* B, int* LDB, int* INFO)
{
    zgesv_(N,NRHS,A,LDA,IPIV,B,LDB,INFO);
}


template <class T> void xtrtrs  (char* UPLO,char* TRANS, char* DIAG, int* N, int* NRHS, const T     * A,int* LDA,T     * B,int* LDB,int* INFO);
template <> void xtrtrs<double> (char* UPLO,char* TRANS, char* DIAG, int* N, int* NRHS, const double* A,int* LDA,double* B,int* LDB,int* INFO)
{
    dtrtrs_(UPLO,TRANS,DIAG,N,NRHS,A,LDA,B,LDB,INFO); //double
}
template <> void xtrtrs<dcmplx> (char* UPLO,char* TRANS, char* DIAG, int* N, int* NRHS, const dcmplx* A,int* LDA,dcmplx* B,int* LDB,int* INFO)
{
    ztrtrs_(UPLO,TRANS,DIAG,N,NRHS,A,LDA,B,LDB,INFO); //complex<double>
}

namespace oml
{
template <class T> typename LapackLinearSolver<T>::VectorT LapackLinearSolver<T>::Solve(const MatrixT& A,const VectorT& b)
{
    assert(A.GetLimits().Row.Low==1);
    assert(A.GetLimits().Col.Low==1);
    int N=A.GetNumRows();
    assert(A.GetNumCols()==N);
    assert(b.size()==N);
    int info=0,nrhs=1;
    Vector<T> x(b);
    Vector<int> ipiv(N);
    MatrixT  A1(A);
    //
    //  Initial call to see how much work space is needed
    //
    xgesv<T>(&N,&nrhs,&A1(1,1),&N,&ipiv(1),&x(1),&N,&info);
    assert(info==0);
    return x;
}

 
template <class T> typename LapackLinearSolver<T>::VectorT LapackLinearSolver<T>::SolveUpperTri(const MatrixT& A,const VectorT& b)
{
    assert(IsUpperTriangular(A));
    return SolveTri(A,b,'U');
}

template <class T> typename LapackLinearSolver<T>::VectorT LapackLinearSolver<T>::SolveLowerTri(const MatrixT& A,const VectorT& b)
{
    assert(IsLowerTriangular(A));
    return SolveTri(A,b,'L');
}

template <class T> typename LapackLinearSolver<T>::VectorT LapackLinearSolver<T>::SolveTri(const MatrixT& A,const VectorT& b,char UL)
{
    assert(A.GetLimits().Row.Low==1);
    assert(A.GetLimits().Col.Low==1);
    int N=A.GetNumRows();
    assert(A.GetNumCols()==N);
    assert(b.size()==N);
    int info=0,nrhs=1;
    char cN('N');
    Vector<T> x(b);
    //
    //  Initial call to see how much work space is needed
    //
    xtrtrs<T>(&UL,&cN,&cN,&N,&nrhs,&A(1,1),&N,&x(1),&N,&info);
    assert(info==0);
    return x;
}

template class LapackLinearSolver<double>;
// template class LapackLinearSolver<dcmplx>;

} //namespace
