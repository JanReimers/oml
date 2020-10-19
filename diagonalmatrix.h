// File: diagonalmatrix.h  make a vector look like a diagonal matrix
#ifndef _DiagonalMatrix_H_
#define _DiagonalMatrix_H_

// Copyright (2020), Jan N. Reimers

#include "oml/iterable.h"
#include "oml/tstream.h"
#include "oml/matrixbase.h"
//#include "oml/matsub.h"
//#include "oml/matrixalg.h"
#include "oml/cow.h"
#include "oml/matindex.h"
//#include "oml/rowcol.h"
#include "oml/vector.h"
//----------------------------------------------------------------------------
//
//  Diagonal matrix class.  For now this is restricted to square matrix shape.
//
template <class T> class DiagonalMatrix
  : public Indexable<T,DiagonalMatrix<T>,Diagonal,Abstract,MatrixShape>
  , public MatrixBase
  , public Iterable<T,DiagonalMatrix<T> >
  , public TStreamableObject<DiagonalMatrix<T> >
{
 public:
  typedef Indexable<T,DiagonalMatrix<T>,Diagonal,Abstract,MatrixShape> IndexableT;
  typedef Ref<T,IndexableT,MatrixShape> RefT;

  explicit DiagonalMatrix();
  explicit DiagonalMatrix(index_t n);
  explicit DiagonalMatrix(index_t nl,index_t nh);
  explicit DiagonalMatrix(const Vector<T>&);
  explicit DiagonalMatrix(const DiagonalMatrix&);
  DiagonalMatrix& operator=(const DiagonalMatrix& dm)
  {
     itsData=dm.itsData;
     return *this;
  }

//  template <class A>                 DiagonalMatrix(const Indexable<T,A,Full,Real,MatrixShape>&);
//  template <class A,Store M, Data D> DiagonalMatrix(const Indexable<T,A,M,D,MatrixShape>&);

//  DMatrix& operator=(const DiagonalMatrix&);
//  template <class A>                 DiagonalMatrix& operator=(const Indexable<T,A,Full,Real,MatrixShape>&);
//  template <class A,Store M, Data D> DiagonalMatrix& operator=(const Indexable<T,A,M,D,MatrixShape>&);

  //DMatrix& operator*=(const DMatrix&);

  std::ostream& Write(std::ostream&) const;
  std::istream& Read (std::istream&)      ;
  friend std::ostream& operator<<(std::ostream& os,const DiagonalMatrix& a)
  {
    return os << static_cast<const TStreamableObject<DiagonalMatrix<T> >& >(a);
  }


        T  operator()(index_t,index_t) const;
        T& operator()(index_t,index_t)      ;

  index_t   size  () const; //Required by iterable.
  MatLimits GetLimits() const;

  const Vector<T>& GetDiagonal() const {return itsData;}
  void SetLimits(index_t n,bool preserve=false)
  {
    MatrixBase::SetLimits(MatLimits(n,n));
    itsData.SetLimits(n,preserve);
  }

  DiagonalMatrix& operator*=(const T& s) {itsData*=s;return *this;}

 private:
  friend class Indexable<T,DiagonalMatrix,Diagonal,Abstract,MatrixShape>;
  friend class Iterable<T,DiagonalMatrix>;

const T* Get() const {return itsData.Get();}
      T* Get()       {return itsData.Get();}



  Vector<T> itsData;   //Copy-On-Write array for the data.
};

//-----------------------------------------------------------------------------
//
//  Macro, which expands to an index checking function call,
//  when DEBUG is on
//
#if DEBUG
  #define CHECK(i,j) GetLimits().CheckIndex(i,j)
#else
  #define CHECK(i,j)
#endif


template <class T> inline T DiagonalMatrix<T>::operator()(index_t i,index_t j) const
{
  static T zero(0);
  CHECK(i,j);
  return i==j ? itsData(i) : zero;
}

template <class T> inline T& DiagonalMatrix<T>::operator()(index_t i,index_t j)
{
  CHECK(i,j);
  assert(i==j);  //Cannot allow write off diagonal elements
  return itsData(i);
}

#undef CHECK

template <class T> inline
DiagonalMatrix<T>::DiagonalMatrix()
  : MatrixBase(MatLimits(0,0))
  , itsData   (0)
  {}

template <class T> inline
DiagonalMatrix<T>::DiagonalMatrix(const Vector<T>& v)
  : MatrixBase(MatLimits(v.size(),v.size()))
  , itsData   (v) //Shallow copy
  {}

template <class T> inline
DiagonalMatrix<T>::DiagonalMatrix(const DiagonalMatrix& dm)
  : MatrixBase(dm.GetLimits())
  , itsData   (dm.itsData) //Shallow copy
  {}



template <class T> DiagonalMatrix<T>& operator*=(DiagonalMatrix<T>& a,const DiagonalMatrix<T>& b)
{
    assert(a.GetLimits()==b.GetLimits());
    a.itsData*=b;
    return a;
}

template <class T> inline index_t DiagonalMatrix<T>::size() const
{
  return GetLimits().size();
}

template <class T> inline MatLimits DiagonalMatrix<T>::GetLimits() const
{
  return MatrixBase::GetLimits();
}

#endif //_DiagonalMatrix_H_
