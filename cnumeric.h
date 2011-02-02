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
