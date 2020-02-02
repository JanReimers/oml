// File: DMatrix.H  Matrix class with direct addressing.
#ifndef _DMatrix_H_
#define _DMatrix_H_

// Copyright (1994-2005), Jan N. Reimers

#include "oml/iterable.h"
#include "oml/tstream.h"
#include "oml/matrixbase.h"
#include "oml/matsub.h"
#include "oml/matrixalg.h"
#include "oml/cow.h"
#include "oml/matindex.h"
#include "oml/rowcol.h"

template <class T> class Array;

//----------------------------------------------------------------------------
//
//  Full matrix class with conventional direct subscripting, i.e. no list of
//  column pointers is maintained.
//
template <class T> class DMatrix
  : public Indexable<T,DMatrix<T>,Full,Real,MatrixShape>
  , public MatrixBase
  , public Iterable<T,DMatrix<T> >
  , public TStreamableObject<DMatrix<T> >
{
 public:
  explicit DMatrix(                );
  explicit DMatrix(index_t r, index_t c);
  explicit DMatrix(subsc_t rl,subsc_t rh,subsc_t cl,subsc_t ch);
  explicit DMatrix(const VecLimits& r,const VecLimits& c);
  explicit DMatrix(const MatLimits&);
  explicit DMatrix(const DMatrix& m);
  template <class A>                 DMatrix(const Indexable<T,A,Full,Real,MatrixShape>&);
  template <class A,Store M, Data D> DMatrix(const Indexable<T,A,M,D,MatrixShape>&);

  DMatrix& operator=(const DMatrix&);
  template <class A>                 DMatrix& operator=(const Indexable<T,A,Full,Real,MatrixShape>&);
  template <class A,Store M, Data D> DMatrix& operator=(const Indexable<T,A,M,D,MatrixShape>&);

  //DMatrix& operator*=(const DMatrix&);

  std::ostream& Write(std::ostream&) const;
  std::istream& Read (std::istream&)      ;
  friend std::ostream& operator<<(std::ostream& os,const DMatrix& a)
  {
    return os << static_cast<const TStreamableObject<DMatrix<T> >& >(a);
  }


  const T& operator()(subsc_t,subsc_t) const;
        T& operator()(subsc_t,subsc_t)      ;

  index_t   size  () const; //Required by iterable.
  MatLimits GetLimits() const;

  void SetLimits(const MatLimits&                           , bool preserve=false);
  void SetLimits(index_t r, index_t c                       , bool preserve=false);
  void SetLimits(subsc_t rl,subsc_t rh,subsc_t cl,subsc_t ch, bool preserve=false);
  void SetLimits(const VecLimits& r,const VecLimits& c      , bool preserve=false);

  void ReIndexRows   (const Array<index_t>& index) {::ReIndexRows   (*this,index);}
  void ReIndexColumns(const Array<index_t>& index) {::ReIndexColumns(*this,index);}
  void SwapRows   (subsc_t i,subsc_t j) {::SwapRows   (*this,i,j);}
  void SwapColumns(subsc_t i,subsc_t j) {::SwapColumns(*this,i,j);}
  DMatrix SubMatrix(const MatLimits& lim) const {DMatrix ret(lim);::SubMatrix(ret,*this);return ret;}

  typedef MatrixRow     <T,DMatrix<T>,Full,Real> RowType;
  typedef MatrixColumn  <T,DMatrix<T>,Full,Real> ColType;
  typedef MatrixDiagonal<T,DMatrix<T>,Full,Real> DiagType;

  RowType  GetRow     (subsc_t row) {return RowType (*this,row);}
  ColType  GetColumn  (subsc_t col) {return ColType (*this,col);}
  DiagType GetDiagonal(           ) {return DiagType(*this    );}

  const RowType  GetRow     (subsc_t row) const {return RowType (*this,row);}
  const ColType  GetColumn  (subsc_t col) const {return ColType (*this,col);}
  const DiagType GetDiagonal(           ) const {return DiagType(*this    );}

  typedef typename Iterable <T,DMatrix>::const_iterator  const_iterator ;
  typedef typename Iterable <T,DMatrix>::iterator iterator;

#if DEBUG
  #define CHECK(i,j) assert(itsLimits.CheckIndex(i,j));
#else
  #define CHECK(i,j)
#endif
//-----------------------------------------------------------------------------
//
//  Allows fast L-value access.  Does COW check during construction.
//
  class Subscriptor
  {
  public:
    Subscriptor(Indexable<T,DMatrix,Full,Real,MatrixShape>& m)
      : itsLimits(m.GetLimits())
      , itsPtr(static_cast<DMatrix&>(m).Get())
      {};

    T& operator()(subsc_t i,subsc_t j) {CHECK(i,j);return itsPtr[itsLimits.Offset(i,j)];}

  private:
    MatLimits itsLimits;
    T*        itsPtr;
  };

#undef CHECK
#if DEBUG
  #define CHECK(i) assert(i>=0&&i<itsSize)
#else
  #define CHECK(i)
#endif
  class ArraySubscriptor
  {
   public:
    ArraySubscriptor(Indexable<T,DMatrix,Full,Real,MatrixShape>& a)
      : itsPtr(static_cast<DMatrix*>(&a)->Get())
      , itsSize(a.size())
      {assert(itsPtr);}
    T& operator[](index_t i) {CHECK(i);return itsPtr[i];}
   private:
    T*      itsPtr;
    index_t itsSize;
  };

#undef CHECK

 private:
  friend class Indexable<T,DMatrix,Full,Real,MatrixShape>;
  friend class Iterable<T,DMatrix>;
  friend class Subscriptor;
  friend class ArraySubscriptor;

  T  operator[](index_t i) const;
  T& operator[](index_t i)      ;

  const T* Get() const; //Required by iterable.
        T* Get()      ; //Required by iterable.
  void  Check () const; //Check internal consistency between limits and cow.

  cow_array<T> itsData;   //Copy-On-Write array for the data.
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


template <class T> inline const T& DMatrix<T>::operator()(subsc_t i,subsc_t j) const
{
  CHECK(i,j);
  return itsData.Get()[GetLimits().Offset(i,j)];
}

template <class T> inline T& DMatrix<T>::operator()(subsc_t i,subsc_t j)
{
  CHECK(i,j);
  return itsData.Get()[GetLimits().Offset(i,j)];
}

#undef CHECK

#if DEBUG
  #define CHECK(i) assert(i>=0 && i<itsData.size())
#else
  #define CHECK(i)
#endif

template <class T> inline  T DMatrix<T>::operator[](index_t i) const
{
  CHECK(i);
  return itsData.Get()[i];
}

template <class T> inline T& DMatrix<T>::operator[](index_t i)
{
  CHECK(i);
  return itsData.Get()[i];
}
#undef CHECK

template <class T> template <class A> inline
DMatrix<T>::DMatrix(const Indexable<T,A,Full,Real,MatrixShape>& m)
  : MatrixBase(m.GetLimits        ())
  , itsData   (GetLimits().size())
  {
    ArrayAssign(*this,m); //Use op[].
  }

template <class T> template <class A,Store M, Data D> inline
DMatrix<T>::DMatrix(const Indexable<T,A,M,D,MatrixShape>& m)
  : MatrixBase(m.GetLimits        ())
  , itsData   (GetLimits().size())
  {
    m.GetLimits();
    MatrixAssign(*this,m); //Use op(i,j).
  }


template <class T> template <class A> inline
DMatrix<T>& DMatrix<T>::operator=(const Indexable<T,A,Full,Real,MatrixShape>& m)
{
  if (size()==0) SetLimits(m.GetLimits());
  ArrayAssign(*this,m); //Use op[].
  return *this;
}

template <class T> template <class A,Store M, Data D> inline
DMatrix<T>& DMatrix<T>::operator=(const Indexable<T,A,M,D,MatrixShape>& m)
{
  if (size()==0) SetLimits(m.GetLimits());
  MatrixAssign(*this,m); //Use op(,).
  return *this;
}

template <class T> DMatrix<T>& operator*=(DMatrix<T>& a,const DMatrix<T>& b)
{
    assert(a.GetLimits().Col==b.GetLimits().Row);
    DMatrix<T> temp=a*b;
    a.SetLimits(temp.GetLimits());
    a=temp;
    return a;
}

template <class T> inline index_t DMatrix<T>::size() const
{
  return GetLimits().size();
}

template <class T> inline const T* DMatrix<T>::Get() const
{
  return itsData.Get();
}

template <class T> inline T* DMatrix<T>::Get()
{
  return itsData.Get();
}

template <class T> inline void DMatrix<T>::SetLimits(const MatLimits& lim,bool preserve)
{
  ::SetLimits(*this,lim,preserve);
}

template <class T> inline void DMatrix<T>::SetLimits(index_t r, index_t c , bool preserve)
{
  SetLimits(MatLimits(r,c),preserve);
}

template <class T> inline void DMatrix<T>::SetLimits(subsc_t rl,subsc_t rh,subsc_t cl,subsc_t ch, bool preserve)
{
  SetLimits(MatLimits(rl,rh,cl,ch),preserve);
}

template <class T> inline void DMatrix<T>::SetLimits(const VecLimits& r,const VecLimits& c , bool preserve)
{
  SetLimits(MatLimits(r,c),preserve);
}

template <class T> inline MatLimits DMatrix<T>::GetLimits() const
{
  return MatrixBase::GetLimits();
}

#endif //_DMatrix_H_
