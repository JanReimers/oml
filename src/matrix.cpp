// File: Matrix.cc  Matrix stored directly in memeory as a 1D array.

// Copyright (1994-2005), Jan N. Reimers

#include "oml/matrix.h"
#include "oml/imp/iterable_io.h"
#include <vector>
#include <iostream>
#include <iomanip>
#include <cassert>

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
  : itsData(0)
  {
    CHECK;
  }

template <class T> Matrix<T>::Matrix(index_t r, index_t c)
  : MatrixBase(MatLimits(r,c))
  , itsData   (GetLimits().size()   )
  {
    CHECK;
  }

template <class T> Matrix<T>::Matrix(index_t rl,index_t rh,index_t cl,index_t ch)
  : MatrixBase(MatLimits(rl,rh,cl,ch))
  , itsData   (GetLimits().size()   )
  {
    CHECK;
  }

template <class T> Matrix<T>::Matrix(const VecLimits& r,const VecLimits& c)
  : MatrixBase(MatLimits(r,c))
  , itsData   (GetLimits().size()   )
  {
    CHECK;
  }

template <class T> Matrix<T>::Matrix(const MatLimits& lim)
  : MatrixBase(lim          )
  , itsData   (lim.size())
  {
    CHECK;
  }

template <class T> Matrix<T>::Matrix(const Matrix<T>& m)
  : MatrixBase(m)
  , itsData   (m.itsData)
  {
//    std::cout << "Matrix<T> copy constructor" << std::endl;
    CHECK;
  }


template <class T> Matrix<T>& Matrix<T>::operator=(const Matrix<T>& m)
{
  MatrixBase::operator=(m);
  itsData=m.itsData;
//  std::cout << "Matrix<T> op=" << std::endl;
  return *this;
}


//-----------------------------------------------------------------------------
//
//  Internal consistency check.
//
template <class T> void Matrix<T>::Check() const
{
  assert(GetLimits().Check());
  assert(GetLimits().size()==static_cast<index_t>(itsData.size()));
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

template <class T> std::istream& Matrix<T>::Read(std::istream& is)
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
//-----------------------------------------------------------------------------
//
//  Changing size and/or limits.
//
template <class T> void Matrix<T>::SetLimits(const MatLimits& theLimits, bool preserve)
{
#ifdef DEBUG
  theLimits.Check();
#endif
  if (GetLimits()!=theLimits)
  {
    if (preserve)
    {
        Matrix<T> dest(theLimits);

      index_t newrowlow=theLimits.Row.Low, newrowhigh=theLimits.Row.High;
      index_t newcollow=theLimits.Col.Low, newcolhigh=theLimits.Col.High;
      index_t rowlow =Max(newrowlow ,GetLimits().Row.Low );   //Limits of old and new
      index_t rowhigh=Min(newrowhigh,GetLimits().Row.High);   //data overlap.
      index_t collow =Max(newcollow ,GetLimits().Col.Low );   //Limits of old and new
      index_t colhigh=Min(newcolhigh,GetLimits().Col.High);   //data overlap.

      Subscriptor sdest (dest);         //Destination subscriptor.

      for (index_t i=rowlow;i<=rowhigh;i++)
        for (index_t j=collow;j<=colhigh;j++)
          sdest(i,j)=(*this)(i,j); //Transfer any overlaping data.

			(*this)=dest;
		}
     else
    {
			(*this)=Matrix<T>(theLimits);
    }
  }
}

template <class T,class A, Store M> void ReIndexRows(Indexable<T,A,M,Real,MatrixShape>& m,const std::vector<index_t>& index)
{
  assert(m.GetLimits().GetNumRows()==static_cast<index_t>(index.size()));

  typename std::vector<index_t>::const_iterator i=index.begin();
  A dest(m.GetLimits());
  typename A::Subscriptor sdest(dest);

  for (int row=m.GetLimits().Row.Low;row<=m.GetLimits().Row.High;row++,i++)
  for (int col=m.GetLimits().Col.Low;col<=m.GetLimits().Col.High;col++)
  sdest(row,col)=m(*i+m.GetLimits().Row.Low,col);

  static_cast<A&>(m)=dest;
}

template <class T,class A, Store M> void ReIndexColumns(Indexable<T,A,M,Real,MatrixShape>& m,const std::vector<index_t>& index)
{
  assert(m.GetLimits().GetNumCols()==static_cast<index_t>(index.size()));

  typename std::vector<index_t>::const_iterator i=index.begin();
  A dest(m.GetLimits());
  typename A::Subscriptor sdest(dest);

  for (int col=m.GetLimits().Col.Low;col<=m.GetLimits().Col.High;col++,i++)
    for (int row=m.GetLimits().Row.Low;row<=m.GetLimits().Row.High;row++)
      sdest(row,col)=m(row,*i+m.GetLimits().Col.Low);

  static_cast<A&>(m)=dest;
}

//-----------------------------------------------------------------------------
//
//  Swapping
//

template <class T> void Matrix<T>::SwapRows(index_t i,index_t j)
{
  Subscriptor s(*this);
  for (index_t c : this->cols())
  {
     T temp=s(i,c);
     s(i,c)=s(j,c);
     s(j,c)=temp;
  }
}

template <class T> void Matrix<T>::SwapColumns(index_t i,index_t j)
{
  Subscriptor s(*this);
  for (index_t r : this->rows())
  {
     T temp=s(r,i);
     s(r,i)=s(r,j);
     s(r,j)=temp;
  }
}

template <class T> Matrix<T> Matrix<T>::SubMatrix(const MatLimits& lim) const
{
    Matrix<T> dest(lim);
	assert(dest.GetLimits().Row.Low >=GetLimits().Row.Low );
	assert(dest.GetLimits().Col.Low >=GetLimits().Col.Low );
	assert(dest.GetLimits().Row.High<=GetLimits().Row.High);
	assert(dest.GetLimits().Col.High<=GetLimits().Col.High);
	Subscriptor s(dest);
	index_t rh=dest.GetLimits().Row.High;
	index_t ch=dest.GetLimits().Col.High;
	for (index_t i=dest.GetLimits().Row.Low;i<=rh;i++)
		for (index_t j=dest.GetLimits().Col.Low;j<=ch;j++)
			s(i,j)=(*this)(i,j);
    return dest;
}


#undef CHECK

