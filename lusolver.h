// File: lusolver.h  System solver using LU decomposition.
#ifndef _lusolver_h_
#define _lusolver_h_

// Copyright (1994-2003), Jan N. Reimers

#include "oml/array.h"

template <class T> class Matrix;
template <class T> class Vector;

/*! \class LUSolver lusolver.h oml/lusolver.h
  \brief Solve a system of equation Ax=B using LU decomposition.
 */
template <class T> class LUSolver
{
 public:
  //! Construct from and immediatly decompose only once.
  LUSolver(Matrix<T>& A);
  //! Solve Ax=V, x is output in V.
  void SolveFor(Vector<T>& V) const;
  //! Solve many equations Ax=B, where x and B are Matrixs, x is output in B.
  void SolveFor(Matrix<T>& B) const;
  bool OK() const {return isOK;}

 private:
  Matrix<T>&     itsA;
  Array<index_t> itsIndex;
  bool           isOK;
};

#endif //_lusolver_h_
