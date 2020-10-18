// File: lubksbm.cc

// Modifications for C++ and oml containers Copyright (1994-2003), Jan N. Reimers

#include "oml/array.h"
#include <cmath>
#include <cassert>
#include <iostream>

#ifndef TYPE
#error "lubksbm.cc TYPE not defined"
#endif

//###########################################################################
//
//  LU backsubstitution of a whole matrix.
//

template <class T> void LUBackSub(const Matrix<T>& a, Matrix<T>& B, const Array<index_t>& Index)
{
  assert(a.GetRowLimits()==a.GetColLimits());
  assert(a.GetColLimits()==B.GetRowLimits());
  assert(a.GetRowLow()==1);
  assert(a.GetColLow()==1);
  assert(B.GetRowLow()==1);
  assert(B.GetColLow()==1);


  typename Matrix<T>::Subscriptor b(B);
  index_t n =a.GetNumRows();

  for (index_t isub=1; isub<=n; isub++)
  {
    int i,ii=0,ip,j;
    T sum;
    for (i=1; i<=n; i++)
    {
      ip=Index[i-1];
      sum=b(ip,isub);
      b(ip,isub)=b(i,isub);
      if (ii)
	for (j=ii; j<=i-1; j++) sum -= a(i,j)*b(j,isub);
      else
	if (sum) ii=i;
      b(i,isub)=sum;
    } //i
    for (i=n;i>=1;i--)
    {
      sum=b(i,isub);
      for (j=i+1; j<=n; j++) sum -= a(i,j)*b(j,isub);
      b(i,isub)=sum/a(i,i);
    } // i
  } // isub
}


typedef TYPE Type;
template void LUBackSub(const Matrix<Type>&, Matrix<Type>&, const Array<index_t>&);
