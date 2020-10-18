// File: gaussj.cc

// Modifications for C++ and oml containers Copyright (1994-2003), Jan N. Reimers

#include "oml/array.h"
#include <cmath>
#include <cassert>
#include <iostream>

#ifndef TYPE
#error "gaussj.cc TYPE not defined"
#endif

//###########################################################################
//
//  Gauss-Jordan elimination for linear system solving.
//  Solves A*X=B for X, answer returned in B, A is destroyed.
//
#define SWAP(a,b) {TYPE temp=(a);(a)=(b);(b)=temp;}


template <class T> void GaussJordanSolver(Matrix<T>& A, Matrix<T>& B)
{
  assert(A.GetRowLimits()==A.GetColLimits());
  assert(B.GetRowLimits()==A.GetColLimits());
  assert(A.GetRowLow()==1);
  assert(A.GetColLow()==1);
  assert(B.GetRowLow()==1);
  assert(B.GetColLow()==1);

  index_t n=A.GetNumRows();
  index_t m=B.GetNumCols();

  Array<index_t> INDXC(n+1),INDXR(n+1),IPIV(n+1);
  index_t icol=0,irow=0;

  Fill(IPIV,0L);

  typename Matrix<T>::Subscriptor a(A);
  typename Matrix<T>::Subscriptor b(B);

  typename Array<index_t>::Subscriptor indxc(INDXC);
  typename Array<index_t>::Subscriptor indxr(INDXR);
  typename Array<index_t>::Subscriptor ipiv (IPIV);

  for (index_t i=1;i<=n;i++)
  {
    T big=0.0;
    for (index_t j=1;j<=n;j++)
      if (ipiv[j] != 1)
	for (index_t k=1;k<=n;k++)
	{
	  if (ipiv[k] == 0)
	  {
	    if (fabs(a(j,k)) >= big)
	    {
	      big=fabs(a(j,k));
	      irow=j;
	      icol=k;
	    } //if
	  } //if
	  else
	    if (ipiv[k] > 1)
	    {
	      std::cerr << "GAUSSJ: Singular Matrix" << std::endl;
	      assert(false);
	    }//if ipiv
	}//k
    ++(ipiv[icol]);
    if (irow != icol)
    {
      for (index_t l=1;l<=n;l++) SWAP(a(irow,l),a(icol,l));
      for (index_t l=1;l<=m;l++) SWAP(b(irow,l),b(icol,l));
    } //if
    indxr[i]=irow;
    indxc[i]=icol;
    if (a(icol,icol) == 0.0)
    {
      std::cerr << "GAUSSJ: Singular Matrix-2" << std::endl;
      assert(false);
    }
    T pivinv=1.0/a(icol,icol);
    a(icol,icol)=1.0;
    for (index_t l=1;l<=n;l++) a(icol,l) *= pivinv;
    for (index_t l=1;l<=m;l++) b(icol,l) *= pivinv;
    for (index_t ll=1;ll<=n;ll++)
    if (ll != icol)
    {
      T dum=a(ll,icol);
      a(ll,icol)=0.0;
      for (index_t l=1;l<=n;l++) a(ll,l) -= a(icol,l)*dum;
      for (index_t l=1;l<=m;l++) b(ll,l) -= b(icol,l)*dum;
    } //if
  }
  for (index_t l=n;l>=1;l--)
  {
    if (indxr[l] != indxc[l])
      for (index_t k=1;k<=n;k++) SWAP(a(k,indxr[l]),a(k,indxc[l]));
  } //for l
}

#undef SWAP

typedef TYPE Type;
template void GaussJordanSolver(Matrix<Type>&, Matrix<Type>&);
