// File: invtri.cc

// Modifications for C++ and oml containers Copyright (1994-2003), Jan N. Reimers

#include "oml/smatrix.h"
#include <cmath>
#include <cassert>

#ifndef TYPE
#error "invtri.cc TYPE not defined"
#endif

//###########################################################################
//
//  Invert an upper triangular matrix.
//
//
template <class T> void InvertTriangular(Matrix<T>& A)
{
  assert(A.GetRowLimits()==A.GetColLimits());
  assert(A.GetRowLow()==1);
  assert(A.GetColLow()==1);

  index_t n=A.GetNumRows();

  typename Matrix<T>::Subscriptor a(A);

  for(subsc_t i=1; i<=n; i++) a(i,i)=1.0/a(i,i);
  for(subsc_t j=n; j>=2; j--)
  {
    for (subsc_t i=j-1; i>=1; i--)
    {
      T temp=0.0;
      for(subsc_t k=i+1; k<=j; k++) temp+=a(i,k)*a(k,j);
      a(i,j)=-a(i,i)*temp;
    }
  }
}

//###########################################################################
//
//  Invert a symmetrically stored upper triangular matrix.
//
//
#define UPPER_ONLY
template <class T> void InvertTriangular(SMatrix<T>& A)
{
  assert(A.GetRowLow()==1);
  assert(A.GetColLow()==1);

  index_t n=A.GetNumRows();

  typename SMatrix<T>::Subscriptor a(A);

  for(subsc_t i=1; i<=n; i++) a(i,i)=1.0/a(i,i);
  for(subsc_t j=n; j>=2; j--)
  {
    for (subsc_t i=j-1; i>=1; i--)
    {
      T temp=0.0;
      for(subsc_t k=i+1; k<=j; k++) temp+=a(i,k)*a(k,j);
      a(i,j)=-a(i,i)*temp;
    }
  }
}
#undef UPPER_ONLY

typedef TYPE Type;
template void InvertTriangular(Matrix<Type>&);
template void InvertTriangular(SMatrix<Type>&);
