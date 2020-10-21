// File: Matrix.cc  Matrix stored directly in memeory as a 1D array.

// Copyright (1994-2005), Jan N. Reimers

#include "oml/diagonalmatrix.h"
#include "oml/array.h"
#include "oml/minmax.h"
#include "oml/iterable_io.h"
#include <iostream>
#include <iomanip>
#include <cassert>

#ifndef TYPE
#error "DiagonalMatrix.cc TYPE was not defined"
#endif

//-----------------------------------------------------------------------------
//
//  IO
//
template <class T> std::ostream& DiagonalMatrix<T>::Write(std::ostream& os) const
{
  assert(os);
  if (!this->Pretty())
  {
    os << GetLimits();
    ::Write(os,*this);
  }
  else
  {
    int prec=os.precision();
    int wid =os.width();
    os << std::setw(0) << GetLimits() << std::endl;
    for (index_t i=GetRowLow();i<=GetRowHigh();i++)
    {
      os << "[ ";
      for (index_t j=GetColLow();j<=GetColHigh();j++)
        os << std::setw(wid) << std::setprecision(prec) << (*this)(i,j) << " ";
      os << "]" << std::endl;
    }
  }
  return os;
}

template <class T> std::istream& DiagonalMatrix<T>::Read(std::istream& is)
{
  assert(is);
  MatLimits lim;
  is >> lim;
  if (size()==0)
    SetLimits(lim.size());
   else
    assert(size()==0);

  ::Read(is,*this);
  assert(is);
  return is;
}


//---------------------------------------------------------------------------------
//
//  Make template instance
//
typedef TYPE Type;
typedef DiagonalMatrix<Type> Mat;
const Store MatStore=Diagonal;

template class DiagonalMatrix<Type>;

//#include "oml/matsub.ci"
