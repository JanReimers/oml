// File: SMatrix.cc  Implementation of a symmetric matrix.

// Copyright (1994-2005), Jan N. Reimers

#include "oml/smatrix.h"
#include "oml/imp/minmax.h"
#include "oml/imp/binio.h"
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
template <class T> SMatrix<T>::SMatrix()
  : itsN   (0)
  , itsData(0)
  {}

template <class T> SMatrix<T>::SMatrix(index_t N)
  : MatrixBase(MatLimits(N,N))
  , itsN      (N)
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
	
template <class T> SMatrix<T>::SMatrix(const Matrix<T>& m, double tol)
  : MatrixBase(m.GetLimits())
  , itsN      (GetLimits().GetNumRows() )
  , itsData   (size(itsN))
  {
    CHECK;
    for (index_t i:this->rows())
        for (index_t j:this->cols(i))
        {
            assert(fabs(m(i,j)-conj(m(j,i)))<=tol);
            (*this)(i,j)=m(i,j);
        }
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
    for (index_t i=GetRowLow();i<=GetRowHigh();i++)
    {
      os << "[ ";
      for (index_t j=GetColLow();j<=GetColHigh();j++)
      {
//        if (i<=j)
	  os << std::setw(wid) << std::setprecision(prec) << (*this)(i,j) << " ";
//	else
//	  for (int n=0;n<=wid;n++) os << " ";
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

//-----------------------------------------------------------------------------
//
//  Changing size and/or limits.
//
template <class T> void SMatrix<T>::SetLimits(const MatLimits& theLimits, bool preserve)
{
#ifdef DEBUG
  theLimits.Check();
#endif
  if (GetLimits()!=theLimits)
  {
    if (preserve)
    {
      SMatrix<T> dest(theLimits);

      index_t newrowlow=theLimits.Row.Low, newrowhigh=theLimits.Row.High;
      //index_t newcollow=theLimits.Col.Low;
      index_t newcolhigh=theLimits.Col.High;
      index_t rowlow =Max(newrowlow ,GetLimits().Row.Low );   //Limits of old and new
      index_t rowhigh=Min(newrowhigh,GetLimits().Row.High);   //data overlap.
//      index_t collow =Max(newcollow ,GetLimits().Col.Low );   //Limits of old and new
      index_t colhigh=Min(newcolhigh,GetLimits().Col.High);   //data overlap.

      Subscriptor sdest (dest);         //Destination subscriptor.

      for (index_t i=rowlow;i<=rowhigh;i++)
        for (index_t j=i;j<=colhigh;j++)
          sdest(i,j)=(*this)(i,j); //Transfer any overlaping data.

			(*this)=dest;
		}
     else
    {
			(*this)=SMatrix<T>(theLimits);
    }
  }
}


template <class T> void SMatrix<T>::ReIndexRows(const std::vector<index_t>& index)
{
  assert(GetLimits().GetNumRows()==index.size());

  typename std::vector<index_t>::const_iterator i=index.begin();
  SMatrix<T> dest(GetLimits());
  Subscriptor sdest(dest);

  for (int row=GetLimits().Row.Low;row<=GetLimits().Row.High;row++,i++)
  for (int col=GetLimits().Col.Low;col<=GetLimits().Col.High;col++)
  sdest(row,col)=(*this)(*i+GetLimits().Row.Low,col);

  *this=dest;
}

template <class T> void SMatrix<T>::ReIndexColumns(const std::vector<index_t>& index)
{
  assert(GetLimits().GetNumCols()==index.size());

  typename std::vector<index_t>::const_iterator i=index.begin();
  SMatrix<T> dest(GetLimits());
  Subscriptor sdest(dest);

  for (int col=GetLimits().Col.Low;col<=GetLimits().Col.High;col++,i++)
    for (int row=GetLimits().Row.Low;row<=GetLimits().Row.High;row++)
      sdest(row,col)=(*this)(row,*i+GetLimits().Col.Low);

  *this=dest;
}
//-----------------------------------------------------------------------------
//
//  Swapping
//

template <class T> void SMatrix<T>::SwapRows(index_t i,index_t j)
{
  Subscriptor s(*this);
  for (index_t c : this->cols())
  {
     T temp=s(i,c);
     s(i,c)=s(j,c);
     s(j,c)=temp;
  }
}

template <class T> void SMatrix<T>::SwapColumns(index_t i,index_t j)
{
  Subscriptor s(*this);
  for (index_t r : this->rows())
  {
     T temp=s(r,i);
     s(r,i)=s(r,j);
     s(r,j)=temp;
  }
}
template <class T> SMatrix<T> SMatrix<T>::SubMatrix(const MatLimits& lim) const
{
    SMatrix<T> dest(lim);
	assert(dest.GetLimits().Row.Low >=GetLimits().Row.Low );
	assert(dest.GetLimits().Col.Low >=GetLimits().Col.Low );
	assert(dest.GetLimits().Row.High<=GetLimits().Row.High);
	assert(dest.GetLimits().Col.High<=GetLimits().Col.High);
	Subscriptor s(dest);
	index_t rh=dest.GetLimits().Row.High;
	index_t ch=dest.GetLimits().Col.High;
	for (index_t i=dest.GetLimits().Row.Low;i<=rh;i++)
		for (index_t j=i;j<=ch;j++)
			s(i,j)=(*this)(i,j);
    return dest;
}

