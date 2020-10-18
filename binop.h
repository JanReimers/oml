// File: binop.h  Glommable Expression Templates.
#ifndef _binop_h_
#define _binop_h_

// Copyright (1994-2003), Jan N. Reimers

#include "oml/mixtypes.h"
#include <cmath>

//-----------------------------------------------------------------
//
//  Concrete binary operations.
//
template <class TA, class TB> class OpAdd
{
 public:
  typedef typename BinaryRetType<TA,TB>::RetType RT;
  static inline RT apply(const TA a, const TB b) {return a+b;}
};

template <class TA, class TB> class OpSub
{
 public:
  typedef typename BinaryRetType<TA,TB>::RetType RT;
  static inline RT apply(const TA a, const TB b) {return a-b;}
};

template <class TA, class TB> class OpMul
{
 public:
  typedef typename BinaryRetType<TA,TB,OpMul>::RetType RT;
   static inline RT apply(const TA a, const TB b) {return a*b;}
};

template <class TA, class TB> class OpDiv
{
 public:
  typedef typename BinaryRetType<TA,TB,OpDiv>::RetType RT;
   static inline RT apply(const TA a, const TB b) {return a/b;}
};

template <class TA, class TB> class OpMod
{
 public:
  typedef typename BinaryRetType<TA,TB>::RetType RT;
  static inline RT apply(const TA a, const TB b) {return a%b;}
};

template <class TA, class TB> class Opfmod
{
 public:
  typedef typename BinaryRetType<TA,TB>::RetType RT;
  static inline RT apply(const TA a, const TB b) {return fmod(a,b);}
};

template <class TA, class TB> class Oppow
{
 public:
  typedef typename BinaryRetType<TA,TB>::RetType RT;
  static inline RT apply(const TA a, const TB b) {return pow(a,b);}
};

template <class TA, class TB> class Opatan2
{
 public:
  typedef typename BinaryRetType<TA,TB>::RetType RT;
  static inline RT apply(const TA a, const TB b) {return atan2(a,b);}
};

template <class TA, class TB> class Ophypot
{
 public:
  typedef typename BinaryRetType<TA,TB>::RetType RT;
  static inline RT apply(const TA a, const TB b) {return hypot(a,b);}
};


template <class TA, class TB> class OpEqual
{
 public:
  typedef bool RT;
  static inline RT apply(const TA a, const TB b) {return a==b;}
};

template <class TA, class TB> class OpAnd
{
 public:
  typedef typename BinaryRetType<TA,TB>::RetType RT;
  static inline RT apply(const TA a, const TB b) {return a&&b;}
};

#include "oml/shape.h"
#include "oml/indext.h"
#include "oml/matlimit.h"
//---------------------------------------------------------------
//
//  Temporary binary operation holder.
//
template <class T, class A, class B, class Op, Shape S> class XprBinOp{};

template <class T, class A, class B, class Op> class XprBinOp<T,A,B,Op,ArrayShape>
{
 public:
   XprBinOp(A a, B b) : itsA(a), itsB(b) {};
  ~XprBinOp() {};

  T       operator[](index_t n) const {return Op::apply(itsA[n],itsB[n]);}
  index_t size      (         ) const {return itsA.size();}
 private:
   A itsA;
   B itsB;
};

template <class T, class A, class B, class Op> class XprBinOp<T,A,B,Op,VectorShape>
{
 public:
   XprBinOp(A a, B b) : itsA(a), itsB(b) {};
  ~XprBinOp() {};

  T         operator[](index_t n) const {return Op::apply(itsA[n],itsB[n]);}
  T         operator()(index_t n) const {return Op::apply(itsA(n),itsB(n));}
  index_t   size      (         ) const {return itsA.size();}
  VecLimits GetLimits (         ) const {return itsA.GetLimits();}
 private:
   A itsA;
   B itsB;
};

template <class T, class A, class B, class Op> class XprBinOp<T,A,B,Op,MatrixShape>
{
 public:
   XprBinOp(A a, B b) : itsA(a), itsB(b) {};
  ~XprBinOp() {};

  T         operator[](index_t n          ) const {return Op::apply(itsA[n],itsB[n]);}
  T         operator()(index_t i,index_t j) const {return Op::apply(itsA(i,j),itsB(i,j));}
  index_t   size      (                   ) const {return itsA.size();}
  MatLimits GetLimits (                   ) const {return itsA.GetLimits();}
 private:
   A itsA;
   B itsB;
};


//---------------------------------------------------------
//
//  Sepecify result of mixing Real and Abstract Data types
//  in binary operators.
//
template <Data D1, Data D2> class BinaryData;
template <Data D> class BinaryData<D,D>
{
 public:
  static const Data RetData=D;
};
template <Data D> class BinaryData<D,Abstract>
{
 public:
  static const Data RetData=Abstract;
};
template <Data D> class BinaryData<Abstract,D>
{
 public:
  static const Data RetData=Abstract;
};
template <> class BinaryData<Abstract,Abstract>
{
 public:
  static const Data RetData=Abstract;
};

#endif //_binop_h_
