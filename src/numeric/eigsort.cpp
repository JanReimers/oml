// File: eigsort.cc

// Copyright (1994-2003), Jan N. Reimers

#include "oml/matrix.h"
#include "oml/vector.h"
#include <cmath>
#include <cassert>

#ifndef TYPE
#error "eigsort.cc TYPE not defined"
#endif

//##########################################################################
//
//  Sorts eigen values in ascending order, vectors are correspondingly
//  rearranged
//

template <class T, class M> void EigenSort(M& A, Vector<T>& EigenValues)
{
  assert(A.GetRowLimits()==A.GetColLimits());
  assert(A.GetRowLow()==1);
  assert(A.GetColLow()==1);

  typename M        ::Subscriptor a(A);
  typename Vector<T>::Subscriptor e(EigenValues);

  index_t n=A.GetNumRows(), k;

  for (index_t i=1;i<n;i++)
  {
    T p=e(k=i);
    for (index_t j=i+1;j<=n;j++) if (e(j) <= p) p=e(k=j);
    if (k != i)
    {
      e(k)=e(i); e(i)=p;
      for (index_t j=1;j<=n;j++) {p=a(j,i); a(j,i)=a(j,k); a(j,k)=p;}
    }
  }
}


typedef TYPE Type;
template void EigenSort(Matrix<Type>&, Vector<Type>&);
