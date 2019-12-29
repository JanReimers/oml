// File: lusolver.cc  System solver using LU decomposition.

// Copyright (1994-2003), Jan N. Reimers

#include "oml/lusolver.h"
#include "oml/vector.h"
#include "oml/numeric.h"

#ifndef TYPE
#error "lusolver.cc TYPE was not defined"
#endif

template <class T> LUSolver<T>::LUSolver(Matrix<T>& A)
  : itsA(A)
  , itsIndex(A.GetNumRows())
  , isOK(true)
  {
    double sign;
    isOK=LUDecomp(itsA,itsIndex,sign);
  }

template <class T> void LUSolver<T>::SolveFor(Vector<T>& V) const
{
  LUBackSub(itsA,V,itsIndex);
}

template <class T> void LUSolver<T>::SolveFor(Matrix<T>& M) const
{
  LUBackSub(itsA,M,itsIndex);
}

typedef TYPE Type;

template class LUSolver<Type>;
