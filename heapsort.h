// File: heapsort.h   Sorting routines.
#ifndef _HeapSort_h_
#define _HeapSort_h_

// Copyright (1994-2003), Jan N. Reimers

#include "oml/indexable.h"

/*! \file heapsort.h
  \brief Heap sorting routines for any OML container.
 */

//-----------------------------------------------------------------------
//
//  Heap sort from NR. You need to supply a comparison operator as 
//  the last template argument.
//
template <class T, class A, Store M, Shape S, class CompOp> 
void HeapSort(Indexable<T,A,M,Real,S>& a)
{
  index_t n=a.size();
  if (n<2) return;
  
  index_t l,j,ir,i;
  T rra;

  typename A::ArraySubscriptor ra(a);  //Make an array subsctriptor, its faster.

  l=(n >> 1)+1;
  ir=n;
  for (;;) {
    if (l > 1)
      rra=ra[--l-1];
    else
    {
      rra=ra[ir-1];
      ra[ir-1]=ra[0];
      if (--ir == 1) {ra[0]=rra;return;}
    }
    i=l;
    j=(l << 1);
    while (j <= ir)
    {
      if (j < ir && CompOp::apply(ra[j-1],ra[j]) ) ++j;
      if (CompOp::apply(rra,ra[j-1])) {ra[i-1]=ra[j-1];j+=(i=j);} else j=ir+1;
    }
    ra[i-1]=rra;
  }
}


//---------------------------------------------------------------------------
//
//  Provides usefull member function for arrays with elements for which
//  < and > operators are defined.
//! Sort any OML container into descending order. Uses heap sort algorithm
template <class T, class A, Store M, Shape S> inline void DescendingHeapSort(Indexable<T,A,M,Real,S>& arr)
{
  HeapSort<T,A,M,S,OpGT<T> >(arr);
}

//! Sort any OML container into ascending order. Uses heap sort algorithm
template <class T, class A, Store M, Shape S> inline void AscendingHeapSort(Indexable<T,A,M,Real,S>& arr)
{
  HeapSort<T,A,M,S,OpLT<T> >(arr);
}

#endif //_HeapSort_h_

