// File: indexsort.h   Sorting routines that create and indexing array.
#ifndef _IndexSort_h_
#define _IndexSort_h_

// Copyright (1994-2003), Jan N. Reimers

#include "oml/imp/indexable.h"
#include <vector>

/*! \file indexsort.h
  \brief Heap sort generated re-indexing arrays.

  Routines for creating re-indexing arrays, so one can sort multple containers
  using one array as the sorting key.

  For example:
  \code
  #include "oml/list.h"
  #include "oml/array.h"

  List<double> X(1000),Y(1000);
  FillRandom(X);
  FillRandom(Y);
  Array<int> index=MakeAscendingIndex(X); //order of X is preserved.
  X.ReIndex(index); //Now X's are sorted.
  Y.ReIndex(index); //Y's are re-arranged exactly as the X's were.

  \endcode
 */

//------------------------------------------------------------------------
//
//  Generate Index Table using heap sort from NR
//
template <class T, class Ind, class CompOp>
std::vector<index_t> HeapIndex(const Ind& arrin)
{
  index_t n=arrin.size();
  std::vector<index_t> index(n);

  index_t l,j,ir,indxt,i;
  T q;

  for (j=0;j<n;j++) index[j]=j;
  if (n<2)  return index;

  l=(n >> 1) + 1;
  ir=n;
  for (;;)
  {
    if (l > 1)
      q=arrin[(indxt=index[--l-1])];
    else
    {
      q=arrin[(indxt=index[ir-1])];
      index[ir-1]=index[0];
      if (--ir == 1) {index[0]=indxt;return index;}
    }
    i=l;
    j=l << 1;
    while (j <= ir)
    {
      if (j < ir && CompOp::apply(arrin[index[j-1]],arrin[index[j]]) ) j++;
      if (CompOp::apply(q,arrin[index[j-1]]))
      {
	index[i-1]=index[j-1];
	j += (i=j);
      }
      else j=ir+1;
    }
    index[i-1]=indxt;
  }
  return index;
}

//-----------------------------------------------------------------------
//
//  Provides usefull member function for arrays with elements for which
//  < and > operators are defined.
//
//! Uses heap sort algorithm to calculate a descending index.
template <class T, class A, Store M, Shape S> inline
std::vector<index_t> MakeDescendingIndex(const ArrayIndexable<T,A,M,S>& arr)
{
  return HeapIndex<T,ArrayIndexable<T,A,M,S>,OpGT<T> >(arr);
}

//! Uses heap sort algorithm to calculate an ascending index.
template <class T, class A, Store M, Data D, Shape S> inline
std::vector<index_t> MakeAscendingIndex(const Indexable<T,A,M,D,S>& arr)
{
  return HeapIndex<T,Indexable<T,A,M,D,S>,OpLT<T> >(arr);
}

//! Verify that an OML containers is sorted.
template <class T, class A, Store M, Data D, Shape S> inline
bool IsAscending(const Indexable<T,A,M,D,S>& a)
{
  bool ret=true;
  for (index_t i=1; i<a.size(); i++) ret=ret&&(a[i-1]<=a[i]);
  return ret;
}

//! Verify that an OML containers is sorted.
template <class T, class A, Store M, Data D, Shape S> inline
bool IsDescending(const Indexable<T,A,M,D,S>& a)
{
  bool ret=true;
  for (index_t i=1; i<a.size(); i++) ret=ret&&(a[i-1]>=a[i]);
  return ret;
}

#endif //_IndexSort_h_

