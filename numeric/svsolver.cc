// File: svsolver.cc  System solver using LU decomposition.

// Copyright (1994-2003), Jan N. Reimers

#include "oml/svsolver.h"
#include "oml/numeric.h"
#include <iostream>
#include <cassert>

#ifndef TYPE
#error "svsolver.cc TYPE was not defined"
#endif

template <class T> SVSolver<T>::SVSolver(Matrix<T>& A)
  : itsU(A)
  , itsV(itsU.GetColLimits(),itsU.GetColLimits())
  , itsW(itsU.GetColLimits())
  {
    SVDecomp(itsU,itsW,itsV);
  }

template <class T> Vector<T> SVSolver<T>::SolveFor(const Vector<T>& B) const
{
  assert(itsU.GetRowLimits()==B.GetLimits());
  Vector<double> X(itsU.GetColLimits());
  SVBackSub(itsU,itsW,itsV,B,X);
  return X;
}

template <class T> Vector<T> SVSolver<T>::SolveFor(const Vector<T>& B, T eps)
{
  assert(itsU.GetRowLimits()==B.GetLimits());
  Vector<T> X(itsU.GetColLimits());

  itsWfixed=itsW;
  double minw=Max(itsWfixed)*eps;
  typedef typename Vector<T>::iterator ITER;
  for(ITER i=itsWfixed.begin();i!=itsWfixed.end();i++) if (*i<minw) *i=0.0;

  SVBackSub(itsU,itsWfixed,itsV,B,X);
  return X;
}

typedef TYPE Type;

//template class SVSolver<Type>;
