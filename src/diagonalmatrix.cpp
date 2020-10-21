// File: Matrix.cc  Matrix stored directly in memeory as a 1D array.

// Copyright (1994-2005), Jan N. Reimers

#include "oml/diagonalmatrix.h"
#include "oml/imp/minmax.h"
#include "oml/imp/iterable_io.h"
#include <iostream>
#include <iomanip>
#include <cassert>

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

template <class T> void DiagonalMatrix<T>::ReIndexRows(const std::vector<index_t>& index)
{
    GetDiagonal().ReIndex(index);
}

template <class T> void DiagonalMatrix<T>::ReIndexColumns(const std::vector<index_t>& index)
{
    GetDiagonal().ReIndex(index);
}

template <class T> void DiagonalMatrix<T>::SwapRows(index_t i,index_t j)
{
     T temp=(*this)(i,i);
     (*this)(i)=(*this)(j,j);
     (*this)(j)=temp;
}

template <class T> void DiagonalMatrix<T>::SwapColumns(index_t i,index_t j)
{
     T temp=(*this)(i,i);
     (*this)(i)=(*this)(j,j);
     (*this)(j)=temp;
}


template <class T> DiagonalMatrix<T> DiagonalMatrix<T>::SubMatrix(const MatLimits& lim) const
{
	assert(lim.Row.Low >=GetLimits().Row.Low );
	assert(lim.Row.High<=GetLimits().Row.High);
	DiagonalMatrix<T> dest(lim);
	for (index_t i:dest.rows())
        dest(i)=(*this)(i,i);
    return dest;
}

