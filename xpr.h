// File: xpr.h  Glommable Expression Templates.
#ifndef _xpr_h_
#define _xpr_h_

// Copyright (1994-2003), Jan N. Reimers

#include "oml/indext.h"
#include "oml/shape.h"
#include "oml/veclimit.h"
#include "oml/matlimit.h"

//----------------------------------------------------
//
//  Hold a reference to terminals
//
template <class T, class R,Shape> class Ref {};

template <class T, class R> class Ref<T,R,ArrayShape>
{
 public:
  Ref(const R& r) : itsRef(r) {};
  T       operator[](index_t n) const {return itsRef[n];}
  index_t size      (         ) const {return itsRef.size();}
 private:
  const R& itsRef;
};

template <class T, class R> class Ref<T,R,VectorShape>
{
 public:
  Ref(const R& r) : itsRef(r) {};
  T         operator[](index_t n) const {return itsRef[n];}
  T         operator()(subsc_t n) const {return itsRef(n);}
  index_t   size      (         ) const {return itsRef.size();}
  VecLimits GetLimits (         ) const {return itsRef.GetLimits();}
 private:
  const R& itsRef;
};

template <class T, class R> class Ref<T,R,MatrixShape>
{
 public:
  Ref(const R& r) : itsRef(r) {};
  T         operator[](index_t n          ) const {return itsRef[n];}
  T         operator()(subsc_t i,subsc_t j) const {return itsRef(i,j);}
  index_t   size      (                   ) const {return itsRef.size();}
  MatLimits GetLimits (                   ) const {return itsRef.GetLimits();}
 private:
  const R& itsRef;
};

//----------------------------------------------------
//
//  Make a scalar or constant look like it is indexalble.
//
template <class T, class A, Shape> class Val {};

template <class T, class A> class Val<T,A,ArrayShape>
{
 public:
  explicit Val(const T& v, const A& a) : itsVal(v), itsA(a), Asize(a.size()) {};
  T       operator[](index_t) const {return itsVal;}
  index_t size      (       ) const {return Asize ;}
 private:
  T        itsVal;
  const A& itsA;
  int Asize; //Kludge to get around a long standing g++ BUG.
};

template <class T, class A> class Val<T,A, VectorShape>
{
 public:
  explicit Val(const T& v, const A& a) : itsVal(v), itsA(a), itsALimits(a.GetLimits()) {};
  T         operator[](index_t) const {return itsVal;}
  T         operator()(subsc_t) const {return itsVal;}
  index_t   size      (       ) const {return itsALimits.size();}
  VecLimits GetLimits (       ) const {return itsALimits;}
 private:
  T        itsVal;
  const A& itsA;
  VecLimits itsALimits; //Kludge to get around a long standing g++ BUG.
};

template <class T, class A> class Val<T,A,MatrixShape>
{
 public:
  explicit Val(const T& v, const A& a) : itsVal(v), itsA(a) {};
  T         operator[](index_t        ) const {return itsVal;}
  T         operator()(subsc_t,subsc_t) const {return itsVal;}
  index_t   size      (               ) const {return itsA.size();}
  MatLimits GetLimits (               ) const {return itsA.GetLimits();}
 private:
  T        itsVal;
  const A& itsA;
};

template <class T, class Derived,Store M,Data D,Shape S> class Indexable;

//---------------------------------------------------
//
//  Hold a temporary expression.
//
template <class T, class Expression,Store M,Data D,Shape S> class Xpr;

template <class T, class Expression,Store M,Data D> class Xpr<T,Expression,M,D,ArrayShape>
: public Indexable<T,Xpr<T,Expression,M,D,ArrayShape>,M,D,ArrayShape>
{
 public:
  Xpr(Expression e) : itsExp(e) {};
  Xpr(const Xpr& x) : itsExp(x.itsExp) {};
  ~Xpr() {};

  T       operator[](index_t n) const {return itsExp[n]       ;}
  index_t size      (         ) const {return itsExp.size();}

 private:
  Expression itsExp;
};

template <class T, class Expression,Store M,Data D> class Xpr<T,Expression,M,D,VectorShape>
: public Indexable<T,Xpr<T,Expression,M,D,VectorShape>,M,D,VectorShape>
{
 public:
  Xpr(Expression e) : itsExp(e) {};
  Xpr(const Xpr& x) : itsExp(x.itsExp) {};
  ~Xpr() {};

  T         operator[](index_t n) const {return itsExp[n];}
  T         operator()(subsc_t i) const {return itsExp(i);}
  index_t   size      (         ) const {return itsExp.size();}
  VecLimits GetLimits (         ) const {return itsExp.GetLimits();}
 private:
  Expression itsExp;
};

template <class T, class Expression,Store M,Data D> class Xpr<T,Expression,M,D,MatrixShape>
: public Indexable<T,Xpr<T,Expression,M,D,MatrixShape>,M,D,MatrixShape>
{
 public:
  Xpr(Expression e) : itsExp(e) {};
  Xpr(const Xpr& x) : itsExp(x.itsExp) {};
  ~Xpr() {};

  T         operator[](index_t n          ) const {return itsExp[n];}
  T         operator()(subsc_t i,subsc_t j) const {return itsExp(i,j);}
  index_t   size      (                   ) const {return itsExp.size();}
  MatLimits GetLimits (                   ) const {return itsExp.GetLimits();}

 private:
  Expression itsExp;
};





#endif //_xpr_h_
