module;
#include <cstddef> //For size_t

export module oml.Xpr;
import oml.VecLimits;
import oml.MatLimits;
import oml.Shape;

export
{
//---------------------------------------------------------------------------
//
//  Primary template for the IndexableBase class which provide index iterators.
//  Use partial specialization for each  container shape.
//
template <class Derived, Shape S> class IndexableBase {};

//---------------------------------------------------------------------------
//
//  Primary template for the indexable class.  Use partial specialization
//  for each  container shape.
//
template <class T, class Derived,Store M,Data D,Shape S> class Indexable;

//---------------------------------------------------------------
//
//  Temporary unary operation holder.
//
template <class T, class TR, class A, Shape S> class XprUnary
{};


template <class T, class TR, class A> class XprUnary<T,TR,A,VectorShape>
{
 public:
  XprUnary(const A& a,TR(*f)(const T&)) : itsA(a), itsF(f) {};
  ~XprUnary() {};

  TR        operator()(index_t n) const {return itsF(itsA(n));}
  size_t    size      (         ) const {return itsA.size();}
  VecLimits GetLimits (         ) const {return itsA.GetLimits();}
 private:
   A itsA;
   TR(*itsF)(const T&) ;
};

//----------------------------------------------------
//
//  Hold a reference to terminals
//
template <class T, class R,Shape> class Ref {};

template <class T, class R> class Ref<T,R,VectorShape>
{
 public:
  Ref(const R& r) : itsRef(r) {};
  T         operator()(index_t n) const {return itsRef(n);}
  size_t    size      (         ) const {return itsRef.size();}
  VecLimits GetLimits (         ) const {return itsRef.GetLimits();}
 private:
  const R& itsRef;
};
template <class Derived, Shape S> class IndexableBase;

//----------------------------------------------------
//
//  Make a scalar or constant look like it is indexalble.
//
template <class T, class A, Shape> class Val {};

template <class T, class A> class Val<T,A, VectorShape>
{
 public:
  explicit Val(const T& v, const A& a) : itsVal(v), itsA(a), itsLimits(a.GetLimits()) {};
  T         operator()(index_t) const {return itsVal;}
  size_t    size      (       ) const {return itsLimits.size();} //Apparently itsA is gone so why store it?
  VecLimits GetLimits (       ) const {return itsLimits;} //Apparently itsA.GetLimits() fails at runtime
 private:
  T        itsVal;
  const A& itsA;
  VecLimits itsLimits;
};


//---------------------------------------------------
//
//  Hold a temporary expression.
//
template <class T, class Expression,Store M,Data D,Shape S> class Xpr;

template <class T, class Expression,Store M,Data D> class Xpr<T,Expression,M,D,VectorShape>
: public Indexable<T,Xpr<T,Expression,M,D,VectorShape>,M,Abstract,VectorShape>
{
 public:
  typedef Indexable<T,Xpr<T,Expression,M,D,VectorShape>,M,Abstract,VectorShape> IndexableT;
  typedef Ref<T,IndexableT,VectorShape> RefT;

  Xpr(Expression e) : itsExp(e) {};
  Xpr(const Xpr& x) : itsExp(x.itsExp) {};
  ~Xpr() {};

  T         operator()(index_t i) const {return itsExp(i);}
  index_t   size      (         ) const {return itsExp.size();}
  VecLimits GetLimits (         ) const {return itsExp.GetLimits();}
 private:
  Expression itsExp;
};

template <class TR, class TA, class TB, class A, class B, Shape S> class XprBinary
{};


template <class TR, class TA, class TB, class A, class B> class XprBinary<TR,TA,TB,A,B,VectorShape>
{
 public:
   XprBinary(const A& a,const B& b,TR(*f)(const TA&,const TB&)) : itsA(a), itsB(b), itsF(f) {};
  ~XprBinary() {};

  TR        operator[](index_t n) const {return itsF(itsA[n],itsB[n]);}
  TR        operator()(index_t n) const {return itsF(itsA(n),itsB(n));}
  size_t    size      (         ) const {return itsA.size();}
  VecLimits GetLimits (         ) const {return itsA.GetLimits();}
 private:
   A itsA;
   B itsB;
   TR(*itsF)(const TA&,const TB&);
};


template <class T, class TR, class A> class XprUnary<T,TR,A,MatrixShape>
{
 public:
  XprUnary(const A& a,TR(*f)(const T&)) : itsA(a), itsF(f) {};
  ~XprUnary() {};

  TR        operator()(index_t i,index_t j) const {return itsF(itsA(i,j));}
  size_t    size      (                   ) const {return itsA.size();}
  MatLimits GetLimits (                   ) const {return itsA.GetLimits();}
 private:
   A itsA;
   TR(*itsF)(const T&);
};


template <class T, class R> class Ref<T,R,MatrixShape>
: public IndexableBase<Ref<T,R,MatrixShape>,MatrixShape> //Get index iterators
{
 public:
  Ref(const R& r) : itsRef(r) {};
  T         operator()(index_t i,index_t j) const {return itsRef(i,j);}
  size_t    size      (                   ) const {return itsRef.size();}
  MatLimits GetLimits (                   ) const {return itsRef.GetLimits();}
 private:
  const R& itsRef;
};


template <class T, class A> class Val<T,A,MatrixShape>
{
 public:
  explicit Val(const T& v, const A& a) : itsVal(v), itsA(a),itsLimits(itsA.GetLimits()) {};
  T         operator()(index_t,index_t) const {return itsVal;}
  index_t   size      (               ) const {return itsLimits.size();} //itsA may be gone
  MatLimits GetLimits (               ) const {return itsLimits;} //Apparently itsA.GetLimits() fails at runtime
 private:
  T         itsVal;
  const A&  itsA;
  MatLimits itsLimits;
};


template <class T, class Expression,Store M,Data D> class Xpr<T,Expression,M,D,MatrixShape>
: public Indexable<T,Xpr<T,Expression,M,D,MatrixShape>,M,Abstract,MatrixShape>
{
 public:
  typedef Indexable<T,Xpr<T,Expression,M,D,MatrixShape>,M,Abstract,MatrixShape> IndexableT;
  typedef Ref<T,IndexableT,MatrixShape> RefT;

  Xpr(Expression e) : itsExp(e) {};
  Xpr(const Xpr& x) : itsExp(x.itsExp) {};
  ~Xpr() {};

  T         operator()(index_t i,index_t j) const {return itsExp(i,j);}
  size_t    size      (                   ) const {return itsExp.size();}
  MatLimits GetLimits (                   ) const {return itsExp.GetLimits();}

 private:
  Expression itsExp;
};

}