// File: Matrix.H  Matrix class with direct column major addressing.
#ifndef _Matrix_H_
#define _Matrix_H_

// Copyright (1994-2005), Jan N. Reimers

#include "oml/imp/matrixbase.h"
#include "oml/imp/arrindex.h"
#include "oml/imp/matindex.h"
#include "oml/imp/tstream.h"
#include "oml/imp/cow.h"
#include "oml/imp/rowcol.h"
#include "oml/imp/minmax.h"
#include <vector>

//----------------------------------------------------------------------------
//
//  Full matrix class with conventional direct subscripting, i.e. no list of
//  column pointers is maintained.
//
template <class T> class Matrix
  : public MatrixBase
  , public ArrayIndexable<T,Matrix<T>,Full     ,MatrixShape>
  , public Indexable     <T,Matrix<T>,Full,Real,MatrixShape>
  , public TStreamableObject<Matrix<T> >
{
 public:
  typedef ArrayIndexable<T,Matrix<T>,Full     ,MatrixShape> IterableT;
  typedef      Indexable<T,Matrix<T>,Full,Real,MatrixShape> IndexableT;
  typedef Ref<T,IndexableT,MatrixShape> RefT;

  explicit Matrix(                );
  explicit Matrix(index_t r, index_t c);
  explicit Matrix(index_t rl,index_t rh,index_t cl,index_t ch);
  explicit Matrix(const VecLimits& r,const VecLimits& c);
  explicit Matrix(const MatLimits&);
  Matrix(const Matrix& m);
  Matrix(Matrix& m) : Matrix<T>(const_cast<const Matrix&>(m)) {};
  template <class A>         Matrix(const ArrayIndexable<T,A,Full,MatrixShape>&);
  template <class A,Store M> Matrix(const Indexable<T,A,M,Abstract,MatrixShape>&);
  template <class A,Store M> Matrix(const Indexable<T,A,M,Real    ,MatrixShape>&);

  Matrix& operator=(const Matrix&);
  template <class A>         Matrix& operator=(const ArrayIndexable<T,A,Full,MatrixShape>&);
  template <class A,Store M> Matrix& operator=(const Indexable<T,A,M,Abstract,MatrixShape>&);
  template <class A,Store M> Matrix& operator=(const Indexable<T,A,M,Real    ,MatrixShape>&);

#ifdef OML_MOVE_OPS
  Matrix(Matrix&& m);
  Matrix& operator=(Matrix&&);
#endif


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
  void RemoveRow   (index_t r);
  void RemoveColumn(index_t r);

  void ReIndexRows   (const std::vector<index_t>& index);
  void ReIndexColumns(const std::vector<index_t>& index);
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

  typedef typename IterableT::const_iterator  const_iterator ;
  typedef typename IterableT::iterator iterator;

  Matrix& operator+=(const ArrayIndexable<T,Matrix,Full,MatrixShape>& b) {return ArrayAdd(*this,b);}
  Matrix& operator-=(const ArrayIndexable<T,Matrix,Full,MatrixShape>& b) {return ArraySub(*this,b);}
  template <class B> Matrix& operator+=(const Indexable<T,B,Full,Abstract,MatrixShape>& b)
  {
      if (size()==0)
      {
        SetLimits(b.GetLimits(),false);
        Fill(*this,T(0));
      }
      return MatrixAdd(*this,b);
  }
  template <class B> Matrix& operator-=(const Indexable<T,B,Full,Abstract,MatrixShape>& b)
  {
      if (size()==0)
      {
        SetLimits(b.GetLimits(),false);
        Fill(*this,T(0));
      }
      return MatrixSub(*this,b);
  }


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
      , itsPtr(static_cast<Matrix&>(m).priv_begin())
      {};

    T& operator()(index_t i,index_t j) {CHECK(i,j);return itsPtr[itsLimits.Offset(i,j)];}

  private:
    MatLimits itsLimits;
    T*        itsPtr;
  };
#undef CHECK

 private:
  friend class Indexable<T,Matrix,Full,Real,MatrixShape>;
  friend class ArrayIndexable <T,Matrix,Full     ,MatrixShape>;
  friend class Subscriptor;

  const T* priv_begin() const {return &*itsData.begin();} //Required by ArrayIndexable.
        T* priv_begin()       {return &*itsData.begin();} //Required by ArrayIndexable.
  void  Check () const; //Check internal consistency between limits and cow.

#ifdef OML_USE_STDVEC
   std::vector<T> itsData;
#else
  cow_array<T> itsData;   //Copy-On-Write array for the data.
#endif
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
  return itsData[GetLimits().Offset(i,j)];
}

template <class T> inline T& Matrix<T>::operator()(index_t i,index_t j)
{
  CHECK(i,j);
  return itsData[GetLimits().Offset(i,j)];
}

#undef CHECK

template <class T> template <class A> inline
Matrix<T>::Matrix(const ArrayIndexable<T,A,Full,MatrixShape>& m)
  : MatrixBase(m.GetLimits        ())
  , itsData   (GetLimits().size())
  {
    ArrayAssign(*this,m); //Use op[].
  }

template <class T> template <class A,Store M> inline
Matrix<T>::Matrix(const Indexable<T,A,M,Abstract,MatrixShape>& m)
  : MatrixBase(m.GetLimits        ())
  , itsData   (GetLimits().size())
  {
    //m.GetLimits();
    MatrixAssign(*this,m); //Use op(i,j).
  }
template <class T> template <class A,Store M> inline
Matrix<T>::Matrix(const Indexable<T,A,M,Real,MatrixShape>& m)
  : MatrixBase(m.GetLimits        ())
  , itsData   (GetLimits().size())
  {
    //m.GetLimits();
    MatrixAssign(*this,m); //Use op(i,j).
  }


template <class T> template <class A> inline
Matrix<T>& Matrix<T>::operator=(const ArrayIndexable<T,A,Full,MatrixShape>& m)
{
  SetLimits(m.GetLimits());
  ArrayAssign(*this,m); //Use op[].
  return *this;
}

template <class T> template <class A,Store M> inline
Matrix<T>& Matrix<T>::operator=(const Indexable<T,A,M,Abstract,MatrixShape>& m)
{
  SetLimits(m.GetLimits());
  MatrixAssign(*this,m); //Use op(,).
  return *this;
}

template <class T> template <class A,Store M> inline
Matrix<T>& Matrix<T>::operator=(const Indexable<T,A,M,Real,MatrixShape>& m)
{
  SetLimits(m.GetLimits());
  MatrixAssign(*this,m); //Use op(,).
  return *this;
}

#ifdef OML_MOVE_OPS

template <class T> inline Matrix<T>::Matrix(Matrix<T>&& m)
  : MatrixBase(m)
  , itsData   (std::move(m.itsData))
  {
//    std::cout << "Matrix<T> move constructor m.itsData.size()=" << m.itsData.size() << std::endl;
  }

template <class T> inline Matrix<T>& Matrix<T>::operator=(Matrix<T>&& m)
{
  MatrixBase::operator=(m);
  itsData=std::move(m.itsData);
//  std::cout << "Matrix<T> move op=" << std::endl;
  return *this;
}
#endif


template <class T, class B,Store M, Data D> Matrix<T>& operator*=(Matrix<T>& a,const Indexable<T,B,M,D,MatrixShape>& b)
{
    assert(a.GetLimits().Col==b.GetLimits().Row);
    Matrix<T> temp=a*b;
    a.SetLimits(temp.GetLimits());
    a=std::move(temp);
    return a;
}

template <class T> inline index_t Matrix<T>::size() const
{
  return GetLimits().size();
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

//------------------------------------------------------------------
//
//  Matrix transpose proxy.
//
template <class T, class Mat> class MatrixTranspose
: public Indexable<T,MatrixTranspose<T,Mat>,Full,Abstract,MatrixShape>
{
 public:
  typedef Indexable<T,MatrixTranspose<T,Mat>,Full,Abstract,MatrixShape> IndexableT;
  typedef Ref<T,IndexableT,MatrixShape> RefT;

  MatrixTranspose(Mat m) : itsMatrix(m) {};
  MatrixTranspose(const MatrixTranspose<T,Mat>& mt) : itsMatrix(mt.itsMatrix) {};

  T operator()(index_t i, index_t j) const {return itsMatrix(j,i);}
  index_t size() const {return GetLimits().size();}
  MatLimits GetLimits() const {return MatLimits(itsMatrix.GetLimits().Col,itsMatrix.GetLimits().Row);}

 private:
  Mat itsMatrix;
};

//-------------------------------------------------------------------------
//
//  Transpose operators, returns a proxy.
//
template <class T, class A, Store M, Data D> inline
MatrixTranspose<T,Ref<T,Indexable<T,A,M,D,MatrixShape>,MatrixShape> >
Transpose(const Indexable<T,A,M,D,MatrixShape>& m)
{
	typedef Ref<T,Indexable<T,A,M,D,MatrixShape>,MatrixShape> ref;
  return MatrixTranspose<T,ref>(ref(m));
}

//! Overloaded operator version of Transpose, returns a proxy.
template <class T, class A, Store M, Data D> inline
MatrixTranspose<T,Ref<T,Indexable<T,A,M,D,MatrixShape>,MatrixShape> >
operator~(const Indexable<T,A,M,D,MatrixShape>& m)
{
  typedef Ref<T,Indexable<T,A,M,D,MatrixShape>,MatrixShape> ref;
  return MatrixTranspose<T,ref>(ref(m));
}


//--------------------------------------------------------------------------
//
//  Matrix algebra functions  M*M, M*V V*M
//
//  All of these operators return a proxy for the operation.  So for example a transpose
//  proxy just stores a Matrix reference and op(i,j) just return ref(j,i). The multiplication
//  proxies are little more complicated, but essentially just return a Row(i)*Col(j) dot
//  product.



//----------------------------------------------------------------------
//
//  Multiplication, Matrix * Matrix Proxy.
//
#include <iostream>
template <class T, class A, class B> class MatrixMMOp
: public Indexable<T,MatrixMMOp<T,A,B>,Full,Abstract,MatrixShape>
{
 public:
  typedef Indexable<T,MatrixMMOp<T,A,B>,Full,Abstract,MatrixShape> IndexableT;
  typedef Ref<T,IndexableT,MatrixShape> RefT;

  MatrixMMOp(const A& a, const B& b)
    : itsA(a)
    , itsB(b)
  {
    //std::cout << "Alim="<< itsA.GetLimits() << " Blim=" << itsB.GetLimits()<< std::endl;
    assert(itsA.GetLimits().Col==itsB.GetLimits().Row);
  };
  MatrixMMOp(const MatrixMMOp& m)
    : itsA(m.itsA)
    , itsB(m.itsB)
  {}
  T operator()(index_t i, index_t j) const
  {
    T ret(0);
    for (index_t k:itsA.cols()) ret+=itsA(i,k)*itsB(k,j);
    return ret;
  }
  index_t size() const {return GetLimits().size();}
  MatLimits GetLimits() const {return MatLimits(itsA.GetLimits().Row,itsB.GetLimits().Col);}

 private:
  const A itsA;
  const B itsB;
};

//----------------------------------------------------------------------
//
//  Multiplication, Matrix * Vector Proxy.
//
template <class T, class A, class V> class MatrixMVOp
: public Indexable<T,MatrixMVOp<T,A,V>,Full,Abstract,VectorShape>
{
 public:
  typedef Indexable<T,MatrixMVOp<T,A,V>,Full,Abstract,VectorShape> IndexableT;
  typedef Ref<T,IndexableT,VectorShape> RefT;

  MatrixMVOp(const A& a, const V& v)
    : itsA(a)
    , itsV(v)
  {
    assert(itsA.GetLimits().Col==itsV.GetLimits());
  };
  MatrixMVOp(const MatrixMVOp& m)
    : itsA(m.itsA)
    , itsV(m.itsV)
    {};
  T operator()(index_t i) const
  {
    T ret(0);
    for (index_t k:itsA.cols()) ret+=itsA(i,k)*itsV(k);
    return ret;
  }
  VecLimits GetLimits() const {return itsA.GetLimits().Row;}
  index_t   size  () const {return GetLimits().size();}

 private:
  const A itsA;
  const V itsV;
};

//----------------------------------------------------------------------
//
//  Multiplication, Vector * Matrix Proxy.
//
template <class T, class V, class B> class MatrixVMOp
: public Indexable<T,MatrixVMOp<T,V,B>,Full,Abstract,VectorShape>
{
 public:
  typedef Indexable<T,MatrixVMOp<T,V,B>,Full,Abstract,VectorShape> IndexableT;
  typedef Ref<T,IndexableT,VectorShape> RefT;

  MatrixVMOp(const V& v, const B& b)
    : itsV(v)
    , itsB(b)
  {
    assert(itsV.GetLimits()==itsB.GetLimits().Row);
  };
  MatrixVMOp(const MatrixVMOp& m)
    : itsV(m.itsV)
    , itsB(m.itsB)
    {};
  T operator()(index_t i) const
  {
    T ret(0);
    for (index_t k: itsB.rows()) ret+=itsV(k)*itsB(k,i);
    return ret;
  }
  VecLimits GetLimits() const {return itsB.GetLimits().Col;}
  index_t   size  () const {return GetLimits().size();}

 private:
  const V itsV;
  const B itsB;
};

//---------------------------------------------------------------------
//
//  Matrix * Matrix, returns a matrix instead of a proxy.
//
template <class TA, class TB, class A, class B, Data DA, Data DB> inline
Matrix<TA> operator*(const Matrix<TA>& a,const Indexable<TB,B,Full,DB,MatrixShape>& b)
{
  typedef typename ReturnType<TA,TB>::RetType TR;
  return Matrix<TA>(MatrixMMOp<TR,typename A::RefT,typename B::RefT>(a,b));
}

//---------------------------------------------------------------------
//
//  Matrix * Matrix, returns a proxy.
//
template <class TA, class TB, class A, class B, Data DA, Data DB> inline
auto operator*(const Indexable<TA,A,Full,DA,MatrixShape>& a,const Indexable<TB,B,Full,DB,MatrixShape>& b)
{
  typedef typename ReturnType<TA,TB>::RetType TR;
  return MatrixMMOp<TR,typename A::RefT,typename B::RefT>(a,b);
}


//--------------------------------------------------------------------------------
//
//  Matrix * Vector, returns a proxy.
//
template <class TA, class TB, class A, class B, Store MA, Store MB, Data DA, Data DB> inline
auto operator*(const Indexable<TA,A,MA,DA,MatrixShape>& a,const Indexable<TB,B,MB,DB,VectorShape>& b)
{
  typedef typename ReturnType<TA,TB>::RetType TR;
  return MatrixMVOp<TR,typename A::RefT,typename B::RefT>(a,b);
}

//------------------------------------------------------------
//
//  Vector * Matrix, returns a proxy.
//
template <class TA, class TB, class A, class B, Store MA, Store MB, Data DA, Data DB> inline
auto operator*(const Indexable<TA,A,MA,DA,VectorShape>& a,const Indexable<TB,B,MB,DB,MatrixShape>& b)
{
  typedef typename ReturnType<TA,TB>::RetType TR;
  return MatrixVMOp<TR,typename A::RefT,typename B::RefT>(a,b);
}

//--------------------------------------------------------------
//
//  Other assorted matrix functions.
//

// Create a unit of Kronecker Delta matrix.
template <class T, class A, Store M> inline
void Unit(Indexable<T,A,M,Real,MatrixShape>& m)
{
  A& a(static_cast<A&>(m));
  Fill(a,T(0.0));
  a.GetDiagonal()=T(1);
}

template <class T, class A, Store M, Data D> inline
A operator~(const Indexable<std::complex<T>,A,M,D,MatrixShape>& m)
{
  return conj(Transpose(m));
}

// Check if matrix is symmetric. This allows for no roundoff errors.
template <class T, class A,Data D> inline
bool IsSymmetric(const Indexable<T,A,Full,D,MatrixShape>& m)
{
  return m==Transpose(m);
}

// Check if matrix is symmetric. This allows for some tolerance.
template <class T, class A,Data D> inline
bool IsSymmetric(const Indexable<T,A,Full,D,MatrixShape>& m,double eps)
{
  return Max(fabs(m-Transpose(m)))<=eps;
}

// Check if matrix is anit-symmetric.
template <class T, class A, Data D> inline
bool IsAntiSymmetric(const Indexable<T,A,Full,D,MatrixShape>& m)
{
  return m==-Transpose(m);
}

// Check if complex matrix  is Hermitian.
template <class T, class A, Data D> inline
bool IsHermitian(const Indexable<std::complex<T>,A,Full,D,MatrixShape>& m)
{
  return m==conj(Transpose(m));
}

// Check if complex matrix  is Hermitian to within a specified precision eps.
template <class T, class A, Data D> inline
bool IsHermitian(const Indexable<std::complex<T>,A,Full,D,MatrixShape>& m, double eps)
{
  return Max(fabs(m-conj(Transpose(m))))<=eps;
}

// Check if complex matrix  is normal to within a specified precision eps.
template <class T, class A, Data D> inline
bool IsNormal(const Indexable<T,A,Full,D,MatrixShape>& m, double eps)
{
  Matrix<T> mdagger=~m;
  return Max(fabs(m*mdagger-mdagger*m))<=eps;
}


// Dummy check for real matrices.
template <class T, class A, Data D> inline
bool IsHermitian(const Indexable<T,A,Full,D,MatrixShape>& m, double eps)
{
  return IsSymmetric(m,eps);
}

// Force a Matrix to be symmetric by averaging off diagonal elements.
template <class T, class A>
double MakeSymmetric(Indexable<T,A,Full,Real,MatrixShape>& m)
{
  A temp=Transpose(m);
  double del=Max(fabs(m-temp));
  static_cast<A&>(m)=(m+temp)/T(2);
  return del;
}

template <class T, class A, Data D> inline
double FrobeniusNorm(const Indexable<T,A,Full,D,MatrixShape>& m)
{
    double fnorm=0.0;
    for (index_t i: m.rows())
        for (index_t j: m.cols())
            fnorm+=real(m(i,j)*conj(m(i,j)));
    return sqrt(fnorm);
}

template <class T, class A, Data D> inline
bool IsUnit(const Indexable<T,A,Full,D,MatrixShape>& m)
{
    Matrix<T> U(m.GetLimits());
    Unit(U);
    return Max(fabs(m-U))==0.0;
}

template <class T, class A, Data D> inline
bool IsUnit(const Indexable<T,A,Full,D,MatrixShape>& m,double eps)
{
    Matrix<T> U(m.GetLimits());
    Unit(U);
    return Max(fabs(m-U))<=eps;
}

template <class T, class A, Data D> inline
bool IsDiagonal(const Indexable<T,A,Full,D,MatrixShape>& m,double eps)
{
    double off(0.0);
    for (index_t i:m.rows())
        for (index_t j:m.cols())
            if (i!=j) off=Max(off,fabs(m(i,j)));
    return off<=eps;
}

template <class T, class A, Data D> inline
bool IsDiagonal(const Indexable<T,A,Full,D,MatrixShape>& m)
{
    return IsDiagonal(m,0.0);
}

template <class T, class A, Data D> inline
bool IsLowerTriangular(const Indexable<T,A,Full,D,MatrixShape>& m)
{
    bool ret=true;
    for (index_t i: m.rows())
        for (index_t j:m.cols(i+1))
        {
            ret = ret && (m(i,j)==0.0);
            if (!ret) break;
        }
    return ret;
}

template <class T, class A, Data D> inline
bool IsLowerTriangular(const Indexable<T,A,Full,D,MatrixShape>& m,double eps)
{
    bool ret=true;
    for (index_t i: m.rows())
        for (index_t j:m.cols(i+1))
        {
            ret = ret && (fabs(m(i,j))<=eps);
            if (!ret) break;
        }
    return ret;
}

template <class T, class A, Data D> inline
bool IsUpperTriangular(const Indexable<T,A,Full,D,MatrixShape>& m)
{
    bool ret=true;
    if (m.GetLimits().GetNumRows()!=0 && m.GetLimits().GetNumCols()!=0)
        for (index_t j: m.cols())
            for (index_t i:m.rows(j+1))
            {
                ret = ret && (m(i,j)==0.0);
                if (!ret) break;
            }
    return ret;
}

template <class T, class A, Data D> inline
bool IsUpperTriangular(const Indexable<T,A,Full,D,MatrixShape>& m,double eps)
{
    bool ret=true;
    if (m.GetLimits().GetNumRows()!=0 && m.GetLimits().GetNumCols()!=0)
        for (index_t j: m.cols())
            for (index_t i:m.rows(j+1))
            {
                ret = ret && (fabs(m(i,j))<=eps);
                if (!ret) break;
            }
    return ret;
}
//----------------------------------------------------------------------
//
//  Vector outer product Proxy.
//
template <class T,class V1,class V2> class MatrixOuterProduct
: public Indexable<T,MatrixOuterProduct<T,V1,V2>,Full,Abstract,MatrixShape>
{
 public:
  typedef Indexable<T,MatrixOuterProduct<T,V1,V2>,Full,Abstract,MatrixShape> IndexableT;
  typedef Ref<T,IndexableT,MatrixShape> RefT;

  MatrixOuterProduct(const V1& a, const V2& b) : itsA(a), itsB(b) {};
  MatrixOuterProduct(const MatrixOuterProduct& op) : itsA(op.itsA), itsB(op.itsB) {};

  T operator()(index_t i, index_t j) const {return itsA(i)*itsB(j);}
  MatLimits GetLimits() const {return MatLimits(itsA.GetLimits(),itsB.GetLimits());}
  index_t size() const {return MatLimits().size();}

 private:
  const V1  itsA;
  const V2  itsB;
};

template <class T,class V> class SymMatrixOuterProduct
: public Indexable<T,SymMatrixOuterProduct<T,V>,Symmetric,Abstract,MatrixShape>
{
 public:
  typedef Indexable<T,SymMatrixOuterProduct<T,V>,Symmetric,Abstract,MatrixShape> IndexableT;
  typedef Ref<T,IndexableT,MatrixShape> RefT;

  SymMatrixOuterProduct(const V& v) : itsV(v) {};
  SymMatrixOuterProduct(const SymMatrixOuterProduct& op) : itsV(op.itsV) {};

  T operator()(index_t i, index_t j) const {return itsV(i)*itsV(j);}
  MatLimits GetLimits() const {return MatLimits(itsV.GetLimits(),itsV.GetLimits());}
  index_t size() const {return MatLimits().size();}

 private:
  const V itsV;
};

//--------------------------------------------------------------------------
//
//  Outer product, returns a proxy.
//
template <class T, class A, class B, Store MA, Store MB, Data DA, Data DB> inline
MatrixOuterProduct<T,Ref<T,Indexable<T,A,MA,DA,VectorShape>,VectorShape>,Ref<T,Indexable<T,B,MB,DB,VectorShape>,VectorShape> >
OuterProduct(const Indexable<T,A,MA,DA,VectorShape>& v1,const Indexable<T,B,MB,DB,VectorShape>& v2)
{
  typedef typename A::RefT refa;
  typedef typename B::RefT refb;
  return MatrixOuterProduct<T,refa,refb>(refa(v1),refb(v2));
}

//template <class T, class A, Store M, Data D> inline
//SymMatrixOuterProduct<T,Ref<T,Indexable<T,A,M,D,VectorShape>,VectorShape> >
//OuterProduct(const Indexable<T,A,M,D,VectorShape>& v)
//{
//  typedef Ref<T,Indexable<T,A,M,D,VectorShape>,VectorShape> ref;
//  return SymMatrixOuterProduct<T,ref>(ref(v));
//}


#endif //_Matrix_H_
