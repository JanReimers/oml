// File: lubksb.cc

// Modifications for C++ and oml containers Copyright (1994-2003), Jan N. Reimers

#include "oml/vector.h"
#include <cmath>
#include <cassert>
#include <iostream>

#ifndef TYPE
#error "lubksb.cc TYPE not defined"
#endif

template <class T> void LUBackSub(const Matrix<T>& a, Vector<T>& B,const std::vector<index_t>& Index)
{
  assert(a.GetRowLimits()==a.GetColLimits());
  assert(B.GetLimits   ()==a.GetRowLimits());
  assert(a.GetRowLow()==1);
  assert(a.GetColLow()==1);
  int i,ii=0,ip,j;
  T sum;

  typename Vector<T>::Subscriptor      b(B);
  index_t n=B.size();

  for (i=1;i<=n;i++)
  {
    ip=Index[i-1];
    sum=b(ip);
    b(ip)=b(i);
    if (ii)
      for (j=ii;j<=i-1;j++) sum -= a(i,j)*b(j);
    else
      if (sum) ii=i;
    b(i)=sum;
  }
  for (i=n;i>=1;i--)
  {
    sum=b(i);
    for (j=i+1;j<=n;j++) sum -= a(i,j)*b(j);
    b(i)=sum/a(i,i);
  }
}

typedef TYPE Type;
template void LUBackSub(const Matrix<Type>&, Vector<Type>&,const std::vector<index_t>&);
