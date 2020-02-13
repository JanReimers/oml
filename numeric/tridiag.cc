// File: tridiag.cc

// Modifications for C++ and oml containers Copyright (1994-2003), Jan N. Reimers

#include "oml/matrix.h"
#include "oml/dmatrix.h"
#include "oml/vector.h"
#include <cmath>
#include <cassert>
#include <iostream>

#ifndef TYPE
#error "tridiag.cc TYPE not defined"
#endif


//###########################################################################
//
//  Housholder reduction of a real symmetric matrix, to tri-diagonal form.
//
//  *this	- Input matrix
//  Diagonal	- Output std::vector of diagonal elements
//  OffDiagonal- Output std::vector of off diagonal elements
//

template <class T, class M> void TriDiagonal(M& A, Vector<T>& Diagonal, Vector<T>& OffDiagonal)
{
  assert(A.GetRowLimits()==A.GetColLimits());
  assert(A.GetRowLow()==1);
  assert(A.GetColLow()==1);

  index_t n=A.GetNumRows();
  T       f,g,hh;

  typename M        ::Subscriptor a (A);
  typename Vector<T>::Subscriptor d (   Diagonal);
  typename Vector<T>::Subscriptor od(OffDiagonal);

  for (subsc_t i=n;i>=2;i--)
  {
    subsc_t l=i-1;
    T h=0.0, scale=0.0;
    if (l > 1)
    {
      for (subsc_t k=1;k<=l;k++) scale += fabs(a(i,k));
      if (scale == 0.0)
      {	od(i)=a(i,l);}
      else
      {
	for (subsc_t k=1;k<=l;k++)
	{
	  a(i,k) /= scale;
	  h += a(i,k)*a(i,k);
	}
	f=a(i,l);
	g = f>0 ? -sqrt(h) : sqrt(h);
	od(i)=scale*g;
	h -= f*g;
	a(i,l)=f-g;
	f=0.0;
	for (subsc_t j=1;j<=l;j++)
	{
	  a(j,i)=a(i,j)/h;
	  g=0.0;
	  for (subsc_t k=1  ;k<=j;k++) g += a(j,k)*a(i,k);
	  for (subsc_t k=j+1;k<=l;k++) g += a(k,j)*a(i,k);
	  od(j)=g/h;
	  f += od(j)*a(i,j);
	}
	hh=f/(h+h);
	for (subsc_t j=1;j<=l;j++)
	{
	  f=a(i,j);
	  od(j)=g=od(j)-hh*f;
	  for (subsc_t k=1;k<=j;k++) a(j,k) -= (f*od(k)+g*a(i,k));
	}
      }
    }
    else od(i)=a(i,l);
    d(i)=h;
  }
  d(1)=0.0;
  od(1)=0.0;
  for (subsc_t i=1;i<=n;i++)
  {
    subsc_t l=i-1;
    if (d(i))
    {
      for (subsc_t j=1;j<=l;j++)
      {
	g=0.0;
	for (subsc_t k=1;k<=l;k++) g += a(i,k)*a(k,j);
	for (subsc_t k=1;k<=l;k++) a(k,j) -= g*a(k,i);
      }
    }
    d(i)=a(i,i);
    a(i,i)=1.0;
    for (subsc_t j=1;j<=l;j++) a(j,i)=a(i,j)=0.0;
  }
}

typedef TYPE Type;
typedef Matrix<Type> M;
typedef DMatrix<Type> DM;
template void TriDiagonal( M&,Vector<Type>&,Vector<Type>&);
template void TriDiagonal(DM&,Vector<Type>&,Vector<Type>&);
