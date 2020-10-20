// File: shellsort.h   Sorting routines.
#ifndef _ShellSort_h_
#define _ShellSort_h_

// Copyright (1994-2003), Jan N. Reimers

#include "oml/indexable.h"

/*! \file shellsort.h
  \brief Shell sorting routines for any OML container.
 */

//---------------------------------------------------------------------
//
//  Shell sort from NR
//
#define ALN2I 1.442695022
#define TINY 1.0e-5

template <class T, class A, Store M, Shape S, class CompOp> 
void ShellSort(Indexable<T,A,M,Real,S>& a)
{
  index_t n=a.size();
  if (n<2) return;
   
  index_t m,lognb2;

  typename A::ArraySubscriptor arr(a); //Make an array subsctriptor, its faster.
  
  lognb2=(index_t)(log((double) n)*ALN2I+TINY);
  m=n;
  for (index_t nn=1;nn<=lognb2;nn++)
  {
    m >>= 1;
    for (index_t j=m;j<n;j++)
    {
      index_t i=j-m;
      T t=arr[j];
      while (i >= 0 && CompOp::apply(arr[i],t) ) {arr[i+m]=arr[i];i -= m;}
      arr[i+m]=t;
    }
  }
}

#undef ALN2I
#undef TINY

//-------------------------------------------------------------------------------
//
//  Provides usefull member function for arrays with elements for which
//  < and > operators are defined.
//
//! Sort any OML container into descending order. Uses shell sort algorithm
template <class T, class A, Store M, Shape S> inline 
void DescendingShellSort(Indexable<T,A,M,Real,S>& arr)
{
  ShellSort<T,A,M,S,OpLT<T> >(arr);
}

//! Sort any OML container into ascending order. Uses shell sort algorithm
template <class T, class A, Store M, Shape S> inline 
void AscendingShellSort(Indexable<T,A,M,Real,S>& arr)
{
  ShellSort<T,A,M,S,OpGT<T> >(arr);
}

#endif //_ShellSort_h_

