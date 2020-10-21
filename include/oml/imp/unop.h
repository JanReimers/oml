// File: unop.h  Glommable Expression Templates.
#ifndef _unop_h_
#define _unop_h_

// Copyright (1994-2005), Jan N. Reimers

#include "oml/imp/shape.h"
#include "oml/imp/index_t.h"
#include "oml/imp/matlimit.h"
#include "oml/imp/mixtypes.h"
#include <cmath>


template <class T> class OpLT
{
 public:
   typedef bool RetType;
   static inline bool apply(const T a, const T b) {return a<b;}
};

template <class T> class OpLE
{
 public:
   typedef bool RetType;
   static inline bool apply(const T a, const T b) {return a<=b;}
};

template <class T> class OpGT
{
 public:
   typedef bool RetType;
   static inline bool apply(const T a, const T b) {return a>b;}
};

template <class T> class OpGE
{
 public:
   typedef bool RetType;
   static inline bool apply(const T a, const T b) {return a>=b;}
};

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

  TR         operator[](index_t n) const {return itsF(itsA[n]);}
  TR         operator()(index_t n) const {return itsF(itsA(n));}
  index_t   size   (         ) const {return itsA.size();}
  VecLimits GetLimits (         ) const {return itsA.GetLimits();}
 private:
   A itsA;
   TR(*itsF)(const T&) ;
};


template <class T, class TR, class A> class XprUnary<T,TR,A,MatrixShape>
{
 public:
  XprUnary(const A& a,TR(*f)(const T&)) : itsA(a), itsF(f) {};
  ~XprUnary() {};

  TR         operator[](index_t n          ) const {return itsF(itsA[n]);}
  TR         operator()(index_t i,index_t j) const {return itsF(itsA(i,j));}
  index_t   size   (                   ) const {return itsA.size();}
  MatLimits GetLimits (                   ) const {return itsA.GetLimits();}
 private:
   A itsA;
   TR(*itsF)(const T&);
};
/*
template <class T, class A, class Op, Shape S> class XprUnaryOp
{};

template <class T, class A, class Op> class XprUnaryOp<T,A,Op,ArrayShape>
{
 public:
  XprUnaryOp(const A& a) : itsA(a) {};
  ~XprUnaryOp() {};

  T       operator[](index_t n) const {return Op::apply(itsA[n]);}
  index_t size   (         ) const {return itsA.size();}
 private:
   A itsA;
};

template <class T, class A, class Op> class XprUnaryOp<T,A,Op,VectorShape>
{
 public:
  XprUnaryOp(const A& a) : itsA(a) {};
  ~XprUnaryOp() {};

  T         operator[](index_t n) const {return Op::apply(itsA[n]);}
  T         operator()(index_t n) const {return Op::apply(itsA(n));}
  index_t   size   (         ) const {return itsA.size();}
  VecLimits GetLimits (         ) const {return itsA.GetLimits();}
 private:
   A itsA;
};

template <class T, class A, class Op> class XprUnaryOp<T,A,Op,MatrixShape>
{
 public:
  XprUnaryOp(const A& a) : itsA(a) {};
  ~XprUnaryOp() {};

  T         operator[](index_t n          ) const {return Op::apply(itsA[n]);}
  T         operator()(index_t i,index_t j) const {return Op::apply(itsA(i,j));}
  index_t   size   (                   ) const {return itsA.size();}
  MatLimits GetLimits (                   ) const {return itsA.GetLimits();}
 private:
   A itsA;
};

*/


#endif //_unop_h_
