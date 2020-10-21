// File: Matrix.cc  Matrix stored as pointers to columns.

// Copyright (1994-2005), Jan N. Reimers

#include "oml/matrix.h"
#include "oml/array.h"
#include "oml/iterable_io.h"
#include "oml/minmax.h"
#include <iostream>
#include <iomanip>
#include <cassert>
#include <cstdlib>

#ifndef TYPE
#error "Matrix.cc TYPE was not defined"
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
template <class T> Matrix<T>::Matrix()
  : itsData   (0)
  , itsColumns(0)
  {
  }

template <class T> Matrix<T>::Matrix(index_t r, index_t c)
  : MatrixBase(MatLimits(r,c))
  , itsData   (GetLimits().size()   )
  , itsColumns(0)
  {
    NewColumnPointers();
    CHECK;
  }

template <class T> Matrix<T>::Matrix(index_t rl,index_t rh,index_t cl,index_t ch)
  : MatrixBase(MatLimits(rl,rh,cl,ch))
  , itsData   (GetLimits().size()   )
  , itsColumns(0)
  {
    NewColumnPointers();
    CHECK;
  }

template <class T> Matrix<T>::Matrix(const VecLimits& r,const VecLimits& c)
  : MatrixBase(MatLimits(r,c))
  , itsData   (GetLimits().size()   )
  , itsColumns(0)
  {
    NewColumnPointers();
    CHECK;
  }

template <class T> Matrix<T>::Matrix(const MatLimits& lim)
  : MatrixBase(MatLimits(lim))
  , itsData   (GetLimits().size()   )
  , itsColumns(0)
  {
    NewColumnPointers();
    CHECK;
  }


template <class T> Matrix<T>::Matrix(const Matrix& m)
  : MatrixBase(m)
  , itsData(m.itsData)
  , itsColumns(m.itsColumns) //First call to Get() will fix these guys up if there is a COW.
  {
    CHECK;
  }

template <class T> Matrix<T>& Matrix<T>::operator=(const Matrix<T>& m)
{
  ReleaseColumns();
  MatrixBase::operator=(m);
  itsData=m.itsData;
  itsColumns=m.itsColumns;
  return *this;
}


//-----------------------------------------------------------------------------
//
//  Internal pointer calculation routines.  These are carfully designed to
//  activat delayed copies of the raw data, if lvalues are returned.
//  When lvalues are requested the column pointers must be updated.
//


template <class T> void Matrix<T>::ReleaseColumns()
{
  if (itsData.GetNumOwners()==1 && itsColumns)
  {
    itsColumns+=GetLimits().Col.Low; //undo fixup.
    //std::cout << "deleteing columns at " << (void*)itsColumns << std::endl;
    delete [] itsColumns;
  }
  itsColumns=0;
}

template <class T> T* Matrix<T>::Get()
{
  if (size()>0)
    if (!CheckColumns()) {itsColumns=0;NewColumnPointers();}  //Force COW check.
  return itsData.Get();
}

template <class T> T** Matrix<T>::GetColumns()
{
  if (size()>0)
    if (!CheckColumns()) {itsColumns=0;NewColumnPointers();} //Force COW check.
  return itsColumns;
}

template <class T> void Matrix<T>::NewColumnPointers()
{
  assert(size()>0);
  index_t num_rows=GetLimits().GetNumRows();
  index_t num_cols=GetLimits().GetNumCols();
  index_t row_low =GetLimits().Row.Low;
  index_t col_low =GetLimits().Col.Low;
  if (itsData.GetNumOwners()>1 || !itsColumns)
  {
    itsColumns=new T*[num_cols];
    // std::cout << "new column pointers at " << (void*)itsColumns << std::endl;
    itsColumns-=GetLimits().Col.Low;
  }
  else
  {
    //    std::cout << "NewColumnPointers() reloading old columns!!" << std::endl;
    exit(-1);
  }
  for (index_t i=0;i<num_cols;i++) itsColumns[i+col_low]=(&itsData.Get()[i*num_rows]-row_low);
}



//-----------------------------------------------------------------------------
//
//  Internal consistency check.
//
template <class T> void Matrix<T>::Check() const
{
  assert(GetLimits().Check());
  assert(GetLimits().size()==itsData.size());
  bool ok=true;
  index_t num_rows=GetNumRows();
  for (index_t i=0;i<GetNumCols();i++)
  {
    ok = ok && itsColumns[i+GetLimits().Col.Low]==(&itsData.Get()[i*num_rows]-GetLimits().Row.Low);
    if (!ok) std::cerr << "col pointer problem i=" << i
      << " colptr=" << (void*)itsColumns[i+GetLimits().Col.Low]
      << " data*=" << (void*)(&itsData.Get()[i*num_rows]-GetLimits().Row.Low) << std::endl;
  }
  assert(ok);
}


//-----------------------------------------------------------------------------
//
//  IO
//
template <class T> std::ostream& Matrix<T>::Write(std::ostream& os) const
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
    const Matrix<T>& s(*this);
    os << std::setw(0) << GetLimits() << std::endl;
    for (index_t i=GetRowLow();i<=GetRowHigh();i++)
    {
      os << "[ ";
      for (index_t j=GetColLow();j<=GetColHigh();j++)
        os  << std::setw(wid) << std::setprecision(prec) << s(i,j) << " ";
      os << "]" << std::endl;
    }
  }
  return os;
}

template <class T> std::istream& Matrix<T>::Read(std::istream& is)
{
  assert(is);
  MatLimits lim;
  is >> lim;
  if (size()==0)
    SetLimits(lim);
   else
    assert(GetLimits()==lim);

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
typedef Matrix<Type> Mat;
const Store MatStore=Full;

template class Matrix<Type>;

#include "oml/matsub.ci"

