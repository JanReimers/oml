// File: invmat.cc  matrix inversion stuff.

// Modifications for C++ and oml containers Copyright (1994-2003), Jan N. Reimers

#include "oml/smatrix.h"
#include "oml/numeric.h"
#include <iostream>

#ifndef TYPE
#error "invmat.cc TYPE not defined"
#endif

//----------------------------------------------------------------------------
//
//  Matrix inversion.
//
template <class T> Matrix<T> InvertSquare(const Matrix<T>& m)
{
  Matrix<T> temp(m);
  Array<index_t> SwapIndex(temp.GetNumRows());
  Matrix<T> inv(temp.GetLimits());
  Unit(inv);
  T sign;
  LUDecomp (temp,SwapIndex,sign);
  LUBackSub(temp,inv,SwapIndex);
  return inv;
}


template <class T> Matrix<T> InvertSymmetric(const Matrix<T>& m)
{
  Matrix<T> ret(m);
  Vector<T> normal=Normalize(ret); //Normalize.
  Cholsky(ret);             //Decompose.
  InvertTriangular(ret);    //Invert upper triangular matrix.
  ret=Matrix<T>(~ret*ret);             //Recompose.
  Normalize(ret,normal);      //Unormalize.
  return ret;
}

#define UPPER_ONLY
template <class T> SMatrix<T> InvertSymmetric(const SMatrix<T>& m)
{
  int n=m.GetNumRows();
  SMatrix<T> ret (m.GetLimits());
  Fill(ret,0.0);

  SMatrix<T> temp(m);
  Vector<T> norm=Normalize(temp); //Normalize.
  Cholsky(temp);             //Decompose.
  InvertTriangular(temp);    //Invert upper triangular matrix.
  typename SMatrix<T>::Subscriptor      sr(ret);

  for (int i=1;i<=n;i++)
    for (int j=i;j<=n;j++)
      for (int k=1;k<=i;k++) sr(i,j)+=temp(k,i)*temp(k,j);

  Normalize(ret,norm);      //Unormalize.
  return ret;
}
#undef _UPPER_ONLY

typedef TYPE Type;
template Matrix   <Type> InvertSquare   (const Matrix   <Type>&);
template Matrix   <Type> InvertSymmetric(const Matrix   <Type>&);
template SMatrix<Type> InvertSymmetric(const SMatrix<Type>&);
