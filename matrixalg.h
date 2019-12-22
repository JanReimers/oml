// File: matrixalg.h  Matrix algebra using expression templates.
#ifndef _matrixalg_h_
#define _matrixalg_h_

// Copyright (1994-2003), Jan N. Reimers

#include "oml/matindex.h"
#include "oml/vector.h"
#include "oml/rowcol.h"

//-----------------------------------------------------------------------------------
//  All of these operators return a proxy for the operation.  So for example a transpose
//  proxy just stores a Matrix reference and op(i,j) just return ref(j,i). The multiplication
//  proxies are little more complicated, but essentially just return a Row(i)*Col(j) dot
//  product.


//------------------------------------------------------------------
//
//  Matrix transpose proxy.
//
template <class T, class Mat> class MatrixTranspose
: public Indexable<T,MatrixTranspose<T,Mat>,Full,Abstract,MatrixShape>
{
 public:
  MatrixTranspose(Mat m) : itsMatrix(m) {};
  MatrixTranspose(const MatrixTranspose<T,Mat>& mt) : itsMatrix(mt.itsMatrix) {};

  T operator()(int i, int j) const {return itsMatrix(j,i);}
  MatLimits GetLimits() const {return MatLimits(itsMatrix.GetLimits().Col,itsMatrix.GetLimits().Row);}

 private:
  Mat itsMatrix;
};

//----------------------------------------------------------------------
//
//  Multiplication, Matrix * Matrix Proxy.
//
template <class T, class A, class B> class MatrixMMOp
: public Indexable<T,MatrixMMOp<T,A,B>,Full,Abstract,MatrixShape>
{
 public:
  MatrixMMOp(const A& a, const B& b)
    : itsA(a)
    , itsB(b)
  {
    assert(itsA.GetLimits().Col==itsB.GetLimits().Row);
  };
  MatrixMMOp(const MatrixMMOp& m)
    : itsA(m.itsA)
    , itsB(m.itsB)
  {}
  T operator()(int i, int j) const
  {
    T ret(0);
    subsc_t ch=itsA.GetLimits().Col.High;
    for (subsc_t k=itsA.GetLimits().Col.Low;k<=ch;k++) ret+=itsA(i,k)*itsB(k,j);
    return ret;
  }
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
  T operator()(int i) const
  {
    T ret(0);
    subsc_t ch=itsA.GetLimits().Col.High;
    for (subsc_t k=itsA.GetLimits().Col.Low;k<=ch;k++) ret+=itsA(i,k)*itsV(k);
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
  T operator()(int i) const
  {
    T ret(0);
    subsc_t rh=itsB.GetLimits().Row.High;
    for (subsc_t k=itsB.GetLimits().Row.Low;k<=rh;k++) ret+=itsV(k)*itsB(k,i);
    return ret;
  }
  VecLimits GetLimits() const {return itsB.GetLimits().Col;}
  index_t   size  () const {return GetLimits()>size();}

 private:
  const V itsV;
  const B itsB;
};

//----------------------------------------------------------------------
//
//  Vector outer product Proxy.
//
template <class T,class V1,class V2> class MatrixOuterProduct
: public Indexable<T,MatrixOuterProduct<T,V1,V2>,Full,Abstract,MatrixShape>
{
 public:
  MatrixOuterProduct(const V1& a, const V2& b) : itsA(a), itsB(b) {};
  MatrixOuterProduct(const MatrixOuterProduct& op) : itsA(op.itsA), itsB(op.itsB) {};

  T operator()(int i, int j) const {return itsA(i)*itsB(j);}
  MatLimits GetLimits() const {return MatLimits(itsA.GetLimits(),itsB.GetLimits());}

 private:
  const V1  itsA;
  const V2  itsB;
};

template <class T,class V> class SymMatrixOuterProduct
: public Indexable<T,SymMatrixOuterProduct<T,V>,Symmetric,Abstract,MatrixShape>
{
 public:
  SymMatrixOuterProduct(const V& v) : itsV(v) {};
  SymMatrixOuterProduct(const SymMatrixOuterProduct& op) : itsV(op.itsV) {};

  T operator()(int i, int j) const {return itsV(i)*itsV(j);}
  MatLimits GetLimits() const {return MatLimits(itsV.GetLimits(),itsV.GetLimits());}

 private:
  const V itsV;
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
//  Outer product, returns a proxy.
//
template <class T, class A, class B, Store MA, Store MB, Data DA, Data DB> inline
MatrixOuterProduct<T,Ref<T,Indexable<T,A,MA,DA,VectorShape>,VectorShape>,Ref<T,Indexable<T,B,MB,DB,VectorShape>,VectorShape> >
OuterProduct(const Indexable<T,A,MA,DA,VectorShape>& v1,const Indexable<T,B,MB,DB,VectorShape>& v2)
{
  typedef Ref<T,Indexable<T,A,MA,DA,VectorShape>,VectorShape> refa;
  typedef Ref<T,Indexable<T,B,MB,DB,VectorShape>,VectorShape> refb;
  return MatrixOuterProduct<T,refa,refb>(refa(v1),refb(v2));
}

template <class T, class A, Store M, Data D> inline
SymMatrixOuterProduct<T,Ref<T,Indexable<T,A,M,D,VectorShape>,VectorShape> >
OuterProduct(const Indexable<T,A,M,D,VectorShape>& v)
{
  typedef Ref<T,Indexable<T,A,M,D,VectorShape>,VectorShape> ref;
  return SymMatrixOuterProduct<T,ref>(ref(v));
}

//---------------------------------------------------------------------
//
//  Matrix * Matrix, returns a proxy.
//
template <class T, class A, class B, Store MA, Store MB, Data DA, Data DB> inline
MatrixMMOp<T,Ref<T,Indexable<T,A,MA,DA,MatrixShape>,MatrixShape>,Ref<T,Indexable<T,B,MB,DB,MatrixShape>,MatrixShape> >
operator*(const Indexable<T,A,MA,DA,MatrixShape>& a,const Indexable<T,B,MB,DB,MatrixShape>& b)
{
  typedef Ref<T,Indexable<T,A,MA,DA,MatrixShape>,MatrixShape> refa;
  typedef Ref<T,Indexable<T,B,MB,DB,MatrixShape>,MatrixShape> refb;
  return MatrixMMOp<T,refa,refb>(refa(a),refb(b));
}

//--------------------------------------------------------------------------------
//
//  Matrix * Vector, returns a proxy.
//
template <class T, class A, class B, Store MA, Store MB, Data DA, Data DB> inline
MatrixMVOp<T,Ref<T,Indexable<T,A,MA,DA,MatrixShape>,MatrixShape>,Ref<T,Indexable<T,B,MB,DB,VectorShape>,VectorShape> >
operator*(const Indexable<T,A,MA,DA,MatrixShape>& a,const Indexable<T,B,MB,DB,VectorShape>& b)
{
  typedef Ref<T,Indexable<T,A,MA,DA,MatrixShape>,MatrixShape> refa;
  typedef Ref<T,Indexable<T,B,MB,DB,VectorShape>,VectorShape> refb;
  return MatrixMVOp<T,refa,refb>(refa(a),refb(b));
}

//------------------------------------------------------------
//
//  Vector * Matrix, returns a proxy.
//
template <class T, class A, class B, Store MA, Store MB, Data DA, Data DB> inline
MatrixVMOp<T,Ref<T,Indexable<T,A,MA,DA,VectorShape>,VectorShape>,Ref<T,Indexable<T,B,MB,DB,MatrixShape>,MatrixShape> >
operator*(const Indexable<T,A,MA,DA,VectorShape>& a,const Indexable<T,B,MB,DB,MatrixShape>& b)
{
  typedef Ref<T,Indexable<T,A,MA,DA,VectorShape>,VectorShape> refa;
  typedef Ref<T,Indexable<T,B,MB,DB,MatrixShape>,MatrixShape> refb;
  return MatrixVMOp<T,refa,refb>(refa(a),refb(b));
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
  a.GetDiagonal()=1;
}

// Check if matrix is symmetric. This allows for no roundoff errors.
template <class T, class A,Data D> inline
bool IsSymmetric(const Indexable<T,A,Full,D,MatrixShape>& m)
{
  return m==Transpose(m);
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

// Force a Matrix to be symmetric by averaging off diagonal elements.
template <class T, class A>
double MakeSymmetric(Indexable<T,A,Full,Real,MatrixShape>& m)
{
  A temp=Transpose(m);
  double del=Max(abs(m-temp));
  static_cast<A&>(m)=(m+temp)/T(2);
  return del;
}

// Normalize from supplied normalization std::vector.
template <class T, class A, class B, Store MB, Data DB> inline
void Normalize(Indexable<T,A,Full,Real,MatrixShape>& m, const Indexable<T,B,MB,DB,VectorShape>& n)
{
  static_cast<A&>(m)=DirectMultiply(m,OuterProduct(n,n));
}

// Normalize from supplied normalization std::vector. Complex version.
template <class T, class A, class B, Store MB, Data DB> inline
void Normalize(Indexable<std::complex<T>,A,Full,Real,MatrixShape>& m, const Indexable<T,B,MB,DB,VectorShape>& n)
{
  static_cast<A&>(m)=DirectMultiply(m,OuterProduct(n,n));
}

template <class T, class A, class B, Store MB, Data DB> inline
void Normalize(Indexable<T,A,Symmetric,Real,MatrixShape>& m, const Indexable<T,B,MB,DB,VectorShape>& n)
{
  static_cast<A&>(m)=DirectMultiply(m,OuterProduct(n));
}

template <class T, class A, class B, Store MB, Data DB> inline
void Normalize(Indexable<std::complex<T>,A,Symmetric,Real,MatrixShape>& m, const Indexable<T,B,MB,DB,VectorShape>& n)
{
  static_cast<A&>(m)=DirectMultiply(m,OuterProduct(n));
}

// Self normalixation using sqrt(GetDiagonal()) for normalization std::vector.
template <class T, class A, Store M> inline Vector<T> Normalize(Indexable<T,A,M,Real,MatrixShape>& m)
{
  Vector<T> ret=T(1)/sqrt(MatrixDiagonal<T,A,M,Real>(m));
  Normalize(m,ret);
  return ret;
}


#endif //_matrixalg_h_
