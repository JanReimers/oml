// File: Numerical.H  Numerical method soubroutines
#ifndef _Numerical_H_
#define _Numerical_H_

// Copyright (1994-2003), Jan N. Reimers

#include "oml/indext.h"
#include "oml/array.h"
#include "oml/vector.h"

template <class T> class Matrix;
template <class T> class SMatrix;

//----------------------------------------------------------------------------
/*! \file numeric.h
  \brief Numerical methods as template functions.
*/
//! Cholsky decomposition of a symmetric matrix, that is \b not stored symmetrically.
template <class T> void Cholsky          (Matrix <T>&                         );
//! Cholsky decomposition of a symmetric matrix, that is stored symmetrically.
template <class T> void Cholsky          (SMatrix<T>&                         );
//! Sort columns of eigenstd::vectors based on the order of eigen values.
template <class T> void EigenSort        (Matrix <T>&, Vector<T>& EigenValues);
//! Solve Ax=B using GJ elimintaion. x is returned in A.
template <class T> void GaussJordanSolver(Matrix <T>& A, Matrix<T>& B          );
//! Invert an upper triangular matrix.
template <class T> void InvertTriangular (Matrix <T>&                         );
//! Invert an upper triangular matrix stored symmetrically.
template <class T> void InvertTriangular (SMatrix<T>&                         );
//! Solve tri-diagonal system TDx=A. x is returned in A.
template <class T> void TriDiagonal      (Matrix <T>& A, Vector<T>& Diagonal, Vector<T>& OffDiagonal);
//! Diagonalize matrix A. Eigen std::vectors returned in A, eigen value returned as a new Vector.
template <class T> Vector<T> Diagonalize(Matrix<T>& A);

//! LU back substitution
template <class T> void LUBackSub(const Matrix<T>& A, Vector<T>& B,const Array<index_t>& Index);
//! LU back substitution
template <class T> void LUBackSub(const Matrix<T>& A, Matrix<T>& B,const Array<index_t>& Index);
//! LU decomposition of A=L*U.
template <class T> bool LUDecomp(Matrix<T>& A, Array<index_t>& Index ,T& d);
//! QL decomposition.  Used for eigen problems.
template <class T> void QLDecomp(Matrix<T>& A, Vector<T>& Diagonal, Vector<T>& OffDiagonal);
//! Single value back substitution.
template <class T> void SVBackSub(const Matrix<T>& U, const Vector<T>& W,const Matrix<T>& V, const Vector<T>& B, Vector<T>& X);
//! Single value decomposition.
template <class T, class M> void SVDecomp (M& A,Vector<T>& W,M& V);

//! Invert a square matrix (not neccesarily symmetric).
template <class T> Matrix <T> InvertSquare   (const Matrix <T>&);
//! Invert a symmtric matrix with full sotrage.
template <class T> Matrix <T> InvertSymmetric(const Matrix <T>&);
//! Invert a symmtric matrix with symmetric sotrage.
template <class T> SMatrix<T> InvertSymmetric(const SMatrix<T>&);

//! Solve symmtric tri diagonal system Tu=r where T is a symmetric tridiagonal matrix. Last element of OffDiagonal is not used.
Array<double> SolveTriDiag(const Array<double> & OffDiagonal,
			   const Array<double> & Diagonal,
			   const Array<double> & r);

//! Solve symmtric tri diagonal system Tu=r where T is a tridiagonal matrix. Last element of OffDiagonal is not used.
Vector<double> SolveTriDiag(const Vector<double> & OffDiagonal,
			    const Vector<double> & Diagonal,
			    const Vector<double> & r);

template <class T> Vector<T> Diagonalize(Matrix<T>& m)
{
  assert(m.GetRowLimits()==m.GetColLimits());
  assert(m.GetRowLimits()==m.GetColLimits());
  assert(IsSymmetric(m));

  Vector<T> EigenValues(m.GetRowLimits());
  Vector<T> OffDiagonal(m.GetRowLimits());

  TriDiagonal(m,EigenValues,OffDiagonal);
  QLDecomp   (m,EigenValues,OffDiagonal);
  EigenSort  (m,EigenValues);
  return EigenValues;
}

#endif // _Numerical_H_
