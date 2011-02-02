// File: SMatrix.cc  Implementation of a symmetric matrix.

// Copyright (1994-2005), Jan N. Reimers

#include "oml/smatrix.h"
#include "oml/array.h"
#include "oml/minmax.h"
#include "oml/binio.h"
#include "oml/iterable_io.h"
#include <iostream>
#include <iomanip>
#include <cassert>

#ifndef TYPE
#error "SMatrix.cc TYPE was not defined"
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
template <class T> SMatrix<T>::SMatrix()
  : itsN   (0)
  , itsData(0)
  {}

template <class T> SMatrix<T>::SMatrix(index_t r, index_t c) 
  : MatrixBase(MatLimits(r,c))
  , itsN      (GetLimits().GetNumRows() )
  , itsData   (size(itsN))
  {
    CHECK;
  }

template <class T> SMatrix<T>::SMatrix(subsc_t rl,subsc_t rh,subsc_t cl,subsc_t ch) 
  : MatrixBase(MatLimits(rl,rh,cl,ch))
  , itsN      (GetLimits().GetNumRows() )
  , itsData   (size(itsN))
  {
    CHECK;
  }

template <class T> SMatrix<T>::SMatrix(const VecLimits& r,const VecLimits& c) 
  : MatrixBase(MatLimits(r,c))
  , itsN      (GetLimits().GetNumRows() )
  , itsData   (size(itsN))
  {
    CHECK;
  }

template <class T> SMatrix<T>::SMatrix(const MatLimits& lim) 
  : MatrixBase(lim              )
  , itsN      (GetLimits().GetNumRows() )
  , itsData   (size(itsN))
  {
    CHECK;
  }

template <class T> SMatrix<T>::SMatrix(const SMatrix<T>& m) 
  : MatrixBase(m)
  , itsN      (m.itsN)
  , itsData   (m.itsData)
  {
    CHECK;
	}


template <class T> SMatrix<T>& SMatrix<T>::operator=(const SMatrix<T>& m)
{
  MatrixBase::operator=(m);
  itsData=m.itsData;
  itsN=m.itsN;
  CHECK;
  return *this;
}

#undef CHECK

//-----------------------------------------------------------------------------
//
//  Internal consistency check.
//
template <class T> void SMatrix<T>::Check() const
{
  assert(GetLimits().Check());
  assert(GetRowLimits()==GetColLimits());
  assert(itsN==GetNumRows());
  assert(size(itsN)==itsData.size());
}


//-----------------------------------------------------------------------------
//
//  IO
//
template <class T> std::ostream& SMatrix<T>::Write(std::ostream& os) const
{
  assert(os);
  if (!this->Pretty())
  {
    os << GetLimits();
    if (this->Binary()) 
      BinaryWrite(itsN,os);
    else
      os << itsN << " ";
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
      {
        if (i<=j) 
	  os << std::setw(wid) << std::setprecision(prec) << (*this)(i,j) << " ";
	else
	  for (int n=0;n<=wid;n++) os << " ";
      }
      os << "]" << std::endl;
    }
  }
  return os;
}

template <class T> std::istream& SMatrix<T>::Read(std::istream& is)
{
  assert(is);
  MatLimits lim;
  is >> lim;
  
  if (this->Binary())
    BinaryRead(itsN,is);
  else
    is >> itsN;
  
  if (size()==0)
    SetLimits(lim);
   else
    assert(GetLimits()==lim);

  ::Read(is,*this);
#if DEBUG
	Check();
#endif
  assert(is);
  return is;
}

//---------------------------------------------------------------------------------
//
//  Make template instance
//
typedef TYPE Type;
typedef SMatrix<Type> Mat;
const Store MatStore=Symmetric;

template class SMatrix<Type>;

#include "oml/matsub.ci"
