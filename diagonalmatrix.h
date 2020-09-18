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
  explicit DiagonalMatrix();
  explicit DiagonalMatrix(index_t n);
  explicit DiagonalMatrix(subsc_t nl,subsc_t nh);
  explicit DiagonalMatrix(const Vector<T>&);

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


        T  operator()(subsc_t,subsc_t) const;
        T& operator()(subsc_t,subsc_t)      ;

  index_t   size  () const; //Required by iterable.
  MatLimits GetLimits() const;

//  void SetLimits(const MatLimits&                           , bool preserve=false);
//  void SetLimits(index_t r, index_t c                       , bool preserve=false);
//  void SetLimits(subsc_t rl,subsc_t rh,subsc_t cl,subsc_t ch, bool preserve=false);
//  void SetLimits(const VecLimits& r,const VecLimits& c      , bool preserve=false);

//  void ReIndexRows   (const Array<index_t>& index) {::ReIndexRows   (*this,index);}
//  void ReIndexColumns(const Array<index_t>& index) {::ReIndexColumns(*this,index);}
//  void SwapRows   (subsc_t i,subsc_t j) {::SwapRows   (*this,i,j);}
//  void SwapColumns(subsc_t i,subsc_t j) {::SwapColumns(*this,i,j);}
//  DMatrix SubMatrix(const MatLimits& lim) const {DMatrix ret(lim);::SubMatrix(ret,*this);return ret;}

//  typedef MatrixRow     <T,DMatrix<T>,Full,Real> RowType;
//  typedef MatrixColumn  <T,DMatrix<T>,Full,Real> ColType;
//  typedef MatrixDiagonal<T,DMatrix<T>,Full,Real> DiagType;

//  RowType  GetRow     (subsc_t row) {return RowType (*this,row);}
//  ColType  GetColumn  (subsc_t col) {return ColType (*this,col);}
//  DiagType GetDiagonal(           ) {return DiagType(*this    );}

//  const RowType  GetRow     (subsc_t row) const {return RowType (*this,row);}
//  const ColType  GetColumn  (subsc_t col) const {return ColType (*this,col);}
//  const DiagType GetDiagonal(           ) const {return DiagType(*this    );}

//  typedef typename Iterable <T,DMatrix>::const_iterator  const_iterator ;
//  typedef typename Iterable <T,DMatrix>::iterator iterator;


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


template <class T> inline T DiagonalMatrix<T>::operator()(subsc_t i,subsc_t j) const
{
  static T zero(0);
  CHECK(i,j);
  return i==j ? itsData(i) : zero;
}

template <class T> inline T& DiagonalMatrix<T>::operator()(subsc_t i,subsc_t j)
{
  CHECK(i,j);
  assert(i==j);  //Cannot allow write off diagonal elements
  return itsData(i);
}

#undef CHECK

template <class T> inline
DiagonalMatrix<T>::DiagonalMatrix(const Vector<T>& v)
  : MatrixBase(MatLimits(v.size(),v.size()))
  , itsData   (v) //Shallow copy
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
