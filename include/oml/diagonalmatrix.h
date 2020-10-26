// File: diagonalmatrix.h  make a vector look like a diagonal matrix
#ifndef _DiagonalMatrix_H_
#define _DiagonalMatrix_H_

// Copyright (2020), Jan N. Reimers

#include "oml/imp/iterable.h"
#include "oml/imp/tstream.h"
#include "oml/imp/matrixbase.h"
#include "oml/imp/matsub.h"
#include "oml/imp/matrixalg.h"
#include "oml/imp/cow.h"
#include "oml/imp/matindex.h"
//#include "oml/rowcol.h"
#include "oml/vector.h"
//----------------------------------------------------------------------------
//
//  Diagonal matrix class.  For now this is restricted to square matrix shape.
//
template <class T> class DiagonalMatrix
    : public Indexable<T,DiagonalMatrix<T>,Diagonal,Real,MatrixShape>
    , public MatrixBase
    , public Iterable<T,DiagonalMatrix<T> >
    , public TStreamableObject<DiagonalMatrix<T> >
{
public:
    typedef Indexable<T,DiagonalMatrix<T>,Diagonal,Real,MatrixShape> IndexableT;
    typedef Ref<T,IndexableT,MatrixShape> RefT;

    explicit DiagonalMatrix();
    explicit DiagonalMatrix(index_t n);
    explicit DiagonalMatrix(const VecLimits&);
    explicit DiagonalMatrix(const MatLimits&);
    explicit DiagonalMatrix(const Vector<T>&);
    DiagonalMatrix(const DiagonalMatrix&);
    DiagonalMatrix& operator=(const DiagonalMatrix& dm)
    {
        itsData=dm.itsData;
        MatrixBase::SetLimits(dm.GetLimits());
        return *this;
    }
    DiagonalMatrix& operator=(const Vector<T>& v)
    {
        itsData=v;
        MatrixBase::SetLimits(MatLimits(v.size(),v.size()));
        return *this;
    }

    template <class A>         DiagonalMatrix(const Indexable<T,A,Diagonal,Real,MatrixShape>&);
    template <class A, Data D> DiagonalMatrix(const Indexable<T,A,Diagonal,D,MatrixShape>&);

    template <class A>         DiagonalMatrix& operator=(const Indexable<T,A,Diagonal,Real,MatrixShape>&);
    template <class A, Data D> DiagonalMatrix& operator=(const Indexable<T,A,Diagonal,D,MatrixShape>&);

#ifdef OML_MOVE_OPS
  DiagonalMatrix(DiagonalMatrix&& m);
  DiagonalMatrix& operator=(DiagonalMatrix&&);
#endif


    std::ostream& Write(std::ostream&) const;
    std::istream& Read (std::istream&)      ;
    friend std::ostream& operator<<(std::ostream& os,const DiagonalMatrix& a)
    {
        return os << static_cast<const TStreamableObject<DiagonalMatrix<T> >& >(a);
    }


    T  operator()(index_t,index_t) const;
    T& operator()(index_t        )      ;

//    T operator[](index_t n) const
//    {
//        return itsData[n];
//    }

    index_t   size  () const; //Required by iterable.
    MatLimits GetLimits() const;

    const Vector<T>& GetDiagonal() const {return itsData;}
          Vector<T>& GetDiagonal()       {return itsData;}

    void SetLimits(const MatLimits& lim,bool preserve=false)
    {
        assert(lim.Row==lim.Col);
        MatrixBase::SetLimits(lim);
        itsData.SetLimits(lim.Row,preserve);
    }

    void SetLimits(index_t n,bool preserve=false)
    {
        SetLimits(MatLimits(n,n),preserve);
    }
    void ReIndexRows   (const std::vector<index_t>& index);
    void ReIndexColumns(const std::vector<index_t>& index);
    void SwapRows   (index_t i,index_t j);
    void SwapColumns(index_t i,index_t j);

    DiagonalMatrix SubMatrix(const MatLimits& lim) const;
    DiagonalMatrix SubMatrix(index_t N) const
    {
        return SubMatrix(MatLimits(N,N));
    }

    DiagonalMatrix& operator*=(const T& s)
    {
        itsData*=s;
        return *this;
    }

//-----------------------------------------------------------------------------
//
//  Allows fast L-value access.  Does COW check during construction.
//
  class Subscriptor
  {
  public:
    Subscriptor(Indexable<T,DiagonalMatrix,Diagonal,Real,MatrixShape>& m)
      : itsLimits(m.GetLimits().Row)
      , itsPtr(static_cast<DiagonalMatrix&>(m).priv_begin())
      {};

    T& operator()(index_t i)
    {
        assert(itsLimits.CheckIndex(i));
        return itsPtr[itsLimits.Offset(i)];
    }

  private:
    VecLimits itsLimits;
    T*        itsPtr;
  };


private:
    friend class Indexable<T,DiagonalMatrix,Diagonal,Real,MatrixShape>;
    friend class Iterable<T,DiagonalMatrix>;

    T  operator[](index_t i) const {return priv_begin()[i];}
    T& operator[](index_t i)       {return priv_begin()[i];}

    const T* priv_begin() const {return itsData.priv_begin();}
          T* priv_begin()       {return itsData.priv_begin();}

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

template <class T> inline T& DiagonalMatrix<T>::operator()(index_t i)
{
    CHECK(i,i);
    return itsData(i);
}

#undef CHECK

template <class T> inline
DiagonalMatrix<T>::DiagonalMatrix()
    : DiagonalMatrix<T>(VecLimits(0))
{}

template <class T> inline
DiagonalMatrix<T>::DiagonalMatrix(index_t N)
    : DiagonalMatrix<T>(VecLimits(N))
{}

template <class T> inline
DiagonalMatrix<T>::DiagonalMatrix(const VecLimits& lim)
    : MatrixBase(MatLimits(lim,lim))
    , itsData   (lim)
{}

template <class T> inline
DiagonalMatrix<T>::DiagonalMatrix(const MatLimits& lim)
    : DiagonalMatrix<T>(lim.Row)
{
    assert(lim.Row==lim.Col);
}

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

template <class T> template <class A>    inline      DiagonalMatrix<T>::
DiagonalMatrix(const Indexable<T,A,Diagonal,Real,MatrixShape>& dm)
    : DiagonalMatrix<T>(dm.GetLimits())
{
    ArrayAssign(*this,dm); //Use op[].
}

template <class T> template <class A, Data D> inline DiagonalMatrix<T>::
DiagonalMatrix(const Indexable<T,A,Diagonal,D,MatrixShape>& dm)
    : DiagonalMatrix<T>(dm.GetLimits())
{
    DiagonalAssign(*this,dm); //Use op(i,j).
}

template <class T> template <class A> inline DiagonalMatrix<T>& DiagonalMatrix<T>::
operator=(const Indexable<T,A,Diagonal,Real,MatrixShape>& dm)
{
    SetLimits(dm.GetLimits());
    ArrayAssign(*this,dm); //Use op[].
    return *this;
}
template <class T> template <class A, Data D> inline DiagonalMatrix<T>& DiagonalMatrix<T>::
operator=(const Indexable<T,A,Diagonal,D,MatrixShape>& dm)
{
    SetLimits(dm.GetLimits());
    DiagonalAssign(*this,dm); //Use op(i,j).
    return *this;
}

#ifdef OML_MOVE_OPS

template <class T> inline DiagonalMatrix<T>::DiagonalMatrix(DiagonalMatrix<T>&& m)
  : MatrixBase(m)
  , itsData   (std::move(m.itsData))
  {
//    std::cout << "DiagonalMatrix<T> move constructor m.itsData.size()=" << m.itsData.size() << std::endl;
  }

template <class T> inline DiagonalMatrix<T>& DiagonalMatrix<T>::operator=(DiagonalMatrix<T>&& m)
{
  MatrixBase::operator=(m);
  itsData=std::move(m.itsData);
//  std::cout << "DiagonalMatrix<T> move op=" << std::endl;
  return *this;
}
#endif


//template <class T> DiagonalMatrix<T>& operator*=(DiagonalMatrix<T>& a,const DiagonalMatrix<T>& b)
//{
//    assert(a.GetLimits()==b.GetLimits());
//    a.itsData*=b;
//    return a;
//}

template <class T> inline index_t DiagonalMatrix<T>::size() const
{
    return itsData.size();
}

template <class T> inline MatLimits DiagonalMatrix<T>::GetLimits() const
{
    return MatrixBase::GetLimits();
}

#endif //_DiagonalMatrix_H_
