// File: svsolver.h  System solver using LU decomposition.
#ifndef _svsolver_h_
#define _svsolver_h_

// Copyright (1994-2003), Jan N. Reimers

#include "oml/vector.h"
#include "oml/matrix.h"

/*! \class SVSolver svsolver.h oml/svsolver.h
  \brief Solve a system of equation Ax=B using single value decomposition.
 */
template <class T> class SVSolver
{
 public:
  //! Construct from and immediatly decompose only once.
   SVSolver(DMatrix<T>& A);
  //! Solve Ax=V, x is output in V.
   Vector<T> SolveFor(const Vector<T>& V) const;
  //! Solve Ax=V, x is output in V. eps is the cutoff on singular values.
   Vector<T> SolveFor(const Vector<T>& V, T eps);

   const Vector<T> GetSingularValues() const {return itsW;}
 private:
  DMatrix<T>&     itsU;
  DMatrix<T>      itsV;
  Vector<T>      itsW;
  Vector<T>      itsWfixed;
};

#endif //_lusolver_h_
