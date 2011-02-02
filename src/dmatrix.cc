// File: DMatrix.cc  Matrix stored directly in memeory as a 1D array.

// Copyright (1994-2005), Jan N. Reimers

#include "oml/dmatrix.h"
#include "oml/array.h"
#include "oml/minmax.h"
#include "oml/iterable_io.h"
#include <iostream>
#include <iomanip>
#include <cassert>

#ifndef TYPE
#error "DMatrix.cc TYPE was not defined"
#endif

#ifdef DEBUG
#define CHECK Check()
#else
#define CHECK
#endif

//-----------------------------------------------------------------------------
//
//  Constructors.
//
template <class T> DMatrix<T>::DMatrix()
  : itsData(0)
  {
    CHECK;
  }

template <class T> DMatrix<T>::DMatrix(index_t r, index_t c) 
  : MatrixBase(MatLimits(r,c))
  , itsData   (GetLimits().size()   )
  {
    CHECK;
  }

template <class T> DMatrix<T>::DMatrix(subsc_t rl,subsc_t rh,subsc_t cl,subsc_t ch) 
  : MatrixBase(MatLimits(rl,rh,cl,ch))
  , itsData   (GetLimits().size()   )
  {
    CHECK;
  }

template <class T> DMatrix<T>::DMatrix(const VecLimits& r,const VecLimits& c) 
  : MatrixBase(MatLimits(r,c))
  , itsData   (GetLimits().size()   )
  {
    CHECK;
  }

template <class T> DMatrix<T>::DMatrix(const MatLimits& lim) 
  : MatrixBase(lim          )
  , itsData   (lim.size())
  {
    CHECK;
  }

template <class T> DMatrix<T>::DMatrix(const DMatrix<T>& m) 
  : MatrixBase(m)
  , itsData   (m.itsData)
  {
    CHECK;
  }

template <class T> DMatrix<T>& DMatrix<T>::operator=(const DMatrix<T>& m)
{
  MatrixBase::operator=(m);
  itsData=m.itsData;
  return *this;
}

//-----------------------------------------------------------------------------
//
//  Internal consistency check.
//
template <class T> void DMatrix<T>::Check() const
{
  assert(GetLimits().Check());
  assert(GetLimits().size()==itsData.size());
}


//-----------------------------------------------------------------------------
//
//  IO
//
template <class T> std::ostream& DMatrix<T>::Write(std::ostream& os) const
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
    for (subsc_t i=GetRowLow();i<=GetRowHigh();i++)
    {
      os << "[ ";
      for (subsc_t j=GetColLow();j<=GetColHigh();j++)
        os << std::setw(wid) << std::setprecision(prec) << (*this)(i,j) << " ";
      os << "]" << std::endl;
    }
  }
  return os;
}

template <class T> std::istream& DMatrix<T>::Read(std::istream& is)
{
  assert(is);
  MatLimits lim;
  is >> lim;
  if (size()==0)
    SetLimits(lim);
   else
    assert(size()==0);

  ::Read(is,*this);
  CHECK;
  assert(is);
  return is;
}


//---------------------------------------------------------------------------------
//
//  Make template instance
//
typedef TYPE Type;
typedef DMatrix<Type> Mat;
const Store MatStore=Full;

template class DMatrix<Type>;

#include "oml/matsub.ci"
