// File: submatrix.h  Proxy class emulating part of an Matrix
#ifndef _submatrix_h_
#define _submatrix_h_

// Copyright (1994-2003), Jan N. Reimers

#include "oml/matrix.h"

template <class T> class RefSubMatrix
  : public Indexable<T,RefSubMatrix<T>,Full,Abstract,MatrixShape>
  , public Iterable<T,RefSubMatrix<T> >
  , public MatrixBase
{
 public:
    RefSubMatrix(Matrix<T>& v, index_t rstart, index_t rstop,index_t cstart, index_t cstop)
      : MatrixBase(MatLimits(0,rstop-rstart,0,cstop-cstart))
      , itsRep(v)
      , itsRowOffset(rstart-v.GetRowLow()) , itsColOffset(cstart-v.GetColLow())
    {
      assert(&itsRep);
      assert(itsRowOffset>=0);
      assert(itsColOffset>=0);
    }

  RefSubMatrix(const RefSubMatrix& sv)
    : MatrixBase(sv.GetLimits())
    , itsRep(sv.itsRep)
    , itsRowOffset(sv.itsRowOffset), itsColOffset(sv.itsColOffset)
   {};

  template <class B> RefSubMatrix& operator=(const Indexable<T,B,Full,Real,MatrixShape>& a)
  {
    assert(GetLimits()==a.GetLimits());
    for (index_t i=0;i<this->GetNumRows();i++)
      for (index_t j=0;j<GetNumCols();j++)
	(*this)(i,j)=a(i,j);
    return *this;
  }

  RefSubMatrix& operator=(const RefSubMatrix& a)
  {
    assert(GetLimits()==a.GetLimits());
    for (index_t i=0;i<GetNumRows();i++)
      for (index_t j=0;j<GetNumCols();j++)
	(*this)(i,j)=a(i,j);
    return *this;
  }

  index_t size() const {return MatrixBase::size();}

  T  operator()(index_t i,index_t j) const
    {
      assert(i>=0);
      assert(i<GetNumRows());
      assert(j>=0);
      assert(j<GetNumCols());
      return const_cast<const Matrix<T>&>(itsRep)(itsRowOffset+i,itsColOffset+j);
    }
  T& operator()(index_t i,index_t j)
    {
      assert(i>=0);
      assert(i<GetNumRows());
      assert(j>=0);
      assert(j<GetNumCols());
      return itsRep(itsRowOffset+i,itsColOffset+j);
    }

  MatLimits GetLimits() const {return MatrixBase::GetLimits();}

  void SetLimits(index_t rstart, index_t rstop,index_t cstart, index_t cstop)
  {
    MatrixBase::SetLimits(MatLimits(rstart,rstop,cstart,cstop));
    itsRowOffset=rstart;
    itsColOffset=cstart;
    assert(itsRowOffset>=0);
    assert(itsColOffset>=0);
    assert(GetNumRows()>=0);
    assert(GetNumCols()>=0);
  }


  Array<T> GetColumn(index_t ic) const
  {
    assert(GetLimits().Col.CheckIndex(ic));
    Array<T> ret(GetNumRows());
    for (index_t ir=0;ir<GetNumRows();ir++) ret[ir]=(*this)(ir,ic);
    return ret;
  }
  Array<T> GetDiagonal() const
  {
    assert(GetLimits().Col==GetLimits().Row);
    Array<T> ret(GetNumRows());
    for (index_t ir=0;ir<GetNumRows();ir++) ret[ir]=(*this)(ir,ir);
    return ret;
  }
#if DEBUG
#define CHECK(i,j) assert(i>=0&&i<itsNumRows); assert(j>=0&&j<itsNumCols)
#else
  #define CHECK(i)
#endif
/*   class MatrixSubscriptor  */
/*   { */
/*    public: */
/*     MatrixSubscriptor(Indexable<T,RefSubMatrix,Full,Real,MatrixShape>& a)  */
/*       : itsPtr(static_cast<RefSubMatrix&>(a).Get()), itsSize(a.size()) {assert(itsPtr);} */
/*     T& operator[](index_t i) {CHECK(i);return itsPtr[i];} */
/*    private: */
/*     T*      itsPtr; */
/*     index_t itsSize; */
/*   }; */
#undef CHECK


 private:
  friend class Iterable <T,RefSubMatrix>;
  friend class Indexable<T,RefSubMatrix,Full,Abstract,MatrixShape>;
  friend class MatrixSubscriptor;

//  T  operator[](index_t i) const {return itsRep[i+itsRowOffset*itsRep.GetLimits().GetNumCols()+itsColOffset];} //Required by Indexable

/*   const T* Get() const {return &itsRep[0]+itsOffset;} //Required by iterable. */
/*         T* Get()       {return &itsRep[0]+itsOffset;} //Required by iterable. */

  Matrix<T>& itsRep;
  index_t   itsRowOffset,itsColOffset;
};

template <class T,class M, class A> inline Array<T>
operator*(const Indexable<T,M,Full,Abstract,MatrixShape>& m, const Indexable<T,A,Full,Real,ArrayShape>& a)
{
  index_t Nr=m.GetLimits().GetNumRows();
  index_t Nc=m.GetLimits().GetNumCols();
  assert(a.size()==Nc);
  Array<T> ret(Nr);
  Fill(ret,0.0);
  for (index_t i=0;i<Nr;i++)
    for (index_t j=0;j<Nc;j++)
      ret[i]+=m(i,j)*a[j];
  return ret;
}
template <class T,class M, class A> inline Array<T>
operator*(const Indexable<T,M,Full,Real,MatrixShape>& m, const Indexable<T,A,Full,Real,ArrayShape>& a)
{
  index_t Nr=m.GetLimits().GetNumRows();
  index_t Nc=m.GetLimits().GetNumCols();
  assert(a.size()==Nc);
  Array<T> ret(Nr);
  Fill(ret,0.0);
  for (index_t i=0;i<Nr;i++)
    for (index_t j=0;j<Nc;j++)
      ret[i]+=m(i,j)*a[j];
  return ret;
}

#endif //_submatrix_h_
