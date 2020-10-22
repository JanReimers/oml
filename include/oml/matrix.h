// File: Matrix.H  Matrix class with direct column major addressing.
#ifndef _Matrix_H_
#define _Matrix_H_

// Copyright (1994-2005), Jan N. Reimers

#include "oml/imp/cow.h"
#include "oml/imp/iterable.h"
#include "oml/imp/tstream.h"
#include "oml/imp/matrixbase.h"
#include "oml/imp/matsub.h"
#include "oml/imp/matindex.h"
#include "oml/imp/rowcol.h"
#include "oml/imp/matrixalg.h"
#include <vector>

//----------------------------------------------------------------------------
//
//  Full matrix class with conventional direct subscripting, i.e. no list of
//  column pointers is maintained.
//
template <class T> class Matrix
  : public Indexable<T,Matrix<T>,Full,Real,MatrixShape>
  , public MatrixBase
  , public Iterable<T,Matrix<T> >
  , public TStreamableObject<Matrix<T> >
{
 public:
  typedef Indexable<T,Matrix<T>,Full,Real,MatrixShape> IndexableT;
  typedef Ref<T,IndexableT,MatrixShape> RefT;

  explicit Matrix(                );
  explicit Matrix(index_t r, index_t c);
  explicit Matrix(index_t rl,index_t rh,index_t cl,index_t ch);
  explicit Matrix(const VecLimits& r,const VecLimits& c);
  explicit Matrix(const MatLimits&);
  explicit Matrix(const Matrix& m);
  template <class A>                 Matrix(const Indexable<T,A,Full,Real,MatrixShape>&);
  template <class A,Store M, Data D> Matrix(const Indexable<T,A,M,D,MatrixShape>&);

  Matrix& operator=(const Matrix&);
  template <class A>                 Matrix& operator=(const Indexable<T,A,Full,Real,MatrixShape>&);
  template <class A,Store M, Data D> Matrix& operator=(const Indexable<T,A,M,D,MatrixShape>&);

  //Matrix& operator*=(const Matrix&);

  std::ostream& Write(std::ostream&) const;
  std::istream& Read (std::istream&)      ;
  friend std::ostream& operator<<(std::ostream& os,const Matrix& a)
  {
    return os << static_cast<const TStreamableObject<Matrix<T> >& >(a);
  }


  const T& operator()(index_t,index_t) const;
        T& operator()(index_t,index_t)      ;

  index_t   size  () const; //Required by iterable.
  MatLimits GetLimits() const;

  void SetLimits(const MatLimits&                           , bool preserve=false);
  void SetLimits(index_t r, index_t c                       , bool preserve=false);
  void SetLimits(index_t rl,index_t rh,index_t cl,index_t ch, bool preserve=false);
  void SetLimits(const VecLimits& r,const VecLimits& c      , bool preserve=false);

  void ReIndexRows   (const std::vector<index_t>& index) {::ReIndexRows   (*this,index);}
  void ReIndexColumns(const std::vector<index_t>& index) {::ReIndexColumns(*this,index);}
  void SwapRows   (index_t i,index_t j);
  void SwapColumns(index_t i,index_t j);
  Matrix SubMatrix(const MatLimits& lim) const;

  typedef MatrixRow     <T,Matrix<T>,Full,Real> RowType;
  typedef MatrixColumn  <T,Matrix<T>,Full,Real> ColType;
  typedef MatrixDiagonal<T,Matrix<T>,Full,Real> DiagType;

  RowType  GetRow     (index_t row) {return RowType (*this,row);}
  ColType  GetColumn  (index_t col) {return ColType (*this,col);}
  DiagType GetDiagonal(           ) {return DiagType(*this    );}

  const RowType  GetRow     (index_t row) const {return RowType (*this,row);}
  const ColType  GetColumn  (index_t col) const {return ColType (*this,col);}
  const DiagType GetDiagonal(           ) const {return DiagType(*this    );}

  typedef typename Iterable <T,Matrix>::const_iterator  const_iterator ;
  typedef typename Iterable <T,Matrix>::iterator iterator;


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
    Subscriptor(Indexable<T,Matrix,Full,Real,MatrixShape>& m)
      : itsLimits(m.GetLimits())
      , itsPtr(static_cast<Matrix&>(m).Get())
      {};

    T& operator()(index_t i,index_t j) {CHECK(i,j);return itsPtr[itsLimits.Offset(i,j)];}

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
    ArraySubscriptor(Indexable<T,Matrix,Full,Real,MatrixShape>& a)
      : itsPtr(static_cast<Matrix*>(&a)->Get())
      , itsSize(a.size())
      {assert(itsPtr);}
    T& operator[](index_t i) {CHECK(i);return itsPtr[i];}
   private:
    T*      itsPtr;
    index_t itsSize;
  };

#undef CHECK

 private:
  friend class Indexable<T,Matrix,Full,Real,MatrixShape>;
  friend class Iterable<T,Matrix>;
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


template <class T> inline const T& Matrix<T>::operator()(index_t i,index_t j) const
{
  CHECK(i,j);
  return itsData.Get()[GetLimits().Offset(i,j)];
}

template <class T> inline T& Matrix<T>::operator()(index_t i,index_t j)
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

template <class T> inline  T Matrix<T>::operator[](index_t i) const
{
  CHECK(i);
  return itsData.Get()[i];
}

template <class T> inline T& Matrix<T>::operator[](index_t i)
{
  CHECK(i);
  return itsData.Get()[i];
}
#undef CHECK

template <class T> template <class A> inline
Matrix<T>::Matrix(const Indexable<T,A,Full,Real,MatrixShape>& m)
  : MatrixBase(m.GetLimits        ())
  , itsData   (GetLimits().size())
  {
    ArrayAssign(*this,m); //Use op[].
  }

template <class T> template <class A,Store M, Data D> inline
Matrix<T>::Matrix(const Indexable<T,A,M,D,MatrixShape>& m)
  : MatrixBase(m.GetLimits        ())
  , itsData   (GetLimits().size())
  {
    //m.GetLimits();
    MatrixAssign(*this,m); //Use op(i,j).
  }


template <class T> template <class A> inline
Matrix<T>& Matrix<T>::operator=(const Indexable<T,A,Full,Real,MatrixShape>& m)
{
  SetLimits(m.GetLimits());
  ArrayAssign(*this,m); //Use op[].
  return *this;
}

template <class T> template <class A,Store M, Data D> inline
Matrix<T>& Matrix<T>::operator=(const Indexable<T,A,M,D,MatrixShape>& m)
{
  SetLimits(m.GetLimits());
  MatrixAssign(*this,m); //Use op(,).
  return *this;
}

template <class T> Matrix<T>& operator*=(Matrix<T>& a,const Matrix<T>& b)
{
    assert(a.GetLimits().Col==b.GetLimits().Row);
    Matrix<T> temp=a*b;
    a.SetLimits(temp.GetLimits());
    a=temp;
    return a;
}

template <class T> inline index_t Matrix<T>::size() const
{
  return GetLimits().size();
}

template <class T> inline const T* Matrix<T>::Get() const
{
  return itsData.Get();
}

template <class T> inline T* Matrix<T>::Get()
{
  return itsData.Get();
}



template <class T> inline void Matrix<T>::SetLimits(index_t r, index_t c , bool preserve)
{
  SetLimits(MatLimits(r,c),preserve);
}

template <class T> inline void Matrix<T>::SetLimits(index_t rl,index_t rh,index_t cl,index_t ch, bool preserve)
{
  SetLimits(MatLimits(rl,rh,cl,ch),preserve);
}

template <class T> inline void Matrix<T>::SetLimits(const VecLimits& r,const VecLimits& c , bool preserve)
{
  SetLimits(MatLimits(r,c),preserve);
}

template <class T> inline MatLimits Matrix<T>::GetLimits() const
{
  return MatrixBase::GetLimits();
}

#endif //_Matrix_H_
