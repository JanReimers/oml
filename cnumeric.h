#include "oml/matrix.h"
#include "oml/smatrix.h"
#include "oml/vector.h"

// Copyright (1994-2003), Jan N. Reimers

/*! \file cnumeric.h
  \brief Numerical methods as template functions for complex matricies.
*/
//! eigen values of symmtric tridiagonal matrix by the rational ql method.
void tqlrat(Vector<double>& d, Vector<double>& e2 ,int ierr);
//! Convert Hermitian matrix into tridiagonal for.
void tql2(Vector<double>& d, Vector<double>& e, Matrix<std::complex<double> >& z, int ierr);
//! Reduce an hermitian matrix to a real symmetric tridiagonal matrix using unitary similarity transformations.
void htridi(Matrix<std::complex<double> >& A,
	    Vector<double>& d ,Vector<double>& e,
	    Vector<double>& e2,Vector<std::complex<double> >& tau);
//! Eigenstd::vectors of a complex hermitian matrix by back transforming those of the corresponding real symmetric tridiagonal matrix determined by  htridi.
void htribk(const Matrix<std::complex<double> >& A,const Vector<std::complex<double> >& tau,
	    Matrix<std::complex<double> >& Z);
//! Get eigenvalues and optionally eigen std::vectors of a complex hermitian matrix.
void ch(Matrix<std::complex<double> >& A, Vector<double>& w ,bool matz,int ierr);

//! Diangaonlize a complex hermitian matrix. Eigen vextors return in A.
template <class T> Vector<T> Diagonalize(Matrix<std::complex<T> >& A, const SMatrix<std::complex<T> >& hm)
{
  assert(hm.GetRowLimits()==hm.GetColLimits());

  Vector<T> EigenValues(hm.GetRowLimits());
  A=hm;
  int err=0;
  ch(A,EigenValues,true,err);
  assert(!err);

  return EigenValues;
}

//! Diangaonlize a complex general matrix. Eigen vextors return in A.
template <class T> Vector<T> Diagonalize(Matrix<std::complex<T> >& A)
{
  assert(A.GetRowLimits()==A.GetColLimits());
  Vector<T> EigenValues(A.GetRowLimits());
  int err=0;
  ch(A,EigenValues,true,err);
  assert(!err);

  return EigenValues;
}

//! Fortran SVD algorithm for complex<double>

extern"C" {
typedef std::complex<double> eType;
//
//      A(MxN) = U(MxM) S(N) V*(NxN)
//      M>=N
void csvd_(eType* a, int* mmax, int* nmax, int* m, int* n, int* p, int* nu, int* nv,
double* s, eType* u, eType* v );

}


using std::cout;
using std::endl;
//! sts::complex<double> version of SVDecomp
template <class TM> void CSVDecomp(TM& A, Vector<double>& s, TM& V)
{
    int M=A.GetNumRows();
    int N=A.GetNumCols();

    bool transpose=N>M;
    if (transpose)
    {
        int it=M;M=N;N=it;
        TM temp=Transpose(A);
        A.SetLimits(0,0);
        A=temp; //TODO make a Transpose member function
        V.SetLimits(N,N);
    }
    assert(N<=M); //Required by fortran routine
    assert(M<=1000);  ////Required by fortran routine
    assert(s.size()==N);
    assert(V.GetNumRows()==N);
    assert(V.GetNumCols()==N);
    TM U(M,M);
    int p=0;
//      A(MxN) = U(MxM) S(N) V*(NxN)
//      M>=N
    csvd_ ( &A(1,1), &M, &N, &M, &N, &p, &N, &N, &s(1), &U(1,1), &V(1,1) );
    if (transpose)
    {
        A.SetLimits(0,0);
        A=conj(V); //Store U in A for return
        V.SetLimits(0,0);
        V=conj(U);
        V.SetLimits(M,N,true); //Reshape
    }
    else
    {
        A=U;
        A.SetLimits(M,N,true); //Reshape
    }
}


