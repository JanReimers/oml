// File: unop.h  Glommable Expression Templates.
#ifndef _unop_h_
#define _unop_h_

// Copyright (1994-2005), Jan N. Reimers

#include "oml/mixtypes.h"
#include <cmath>

#if defined(__CYGWIN32__) || (__GNUC__==2)
extern "C" {double pow(double,double);}
inline double pow10(double x) {return pow(10.0,x);}
#endif

namespace std {template <class T> class complex;}


//-----------------------------------------------------------------
//
//  Concrete unary operations.
//
template <class T> class OpPlus
{
 public:
   typedef  T RetType;

   static inline T apply(const T a) {return +a;}
};

template <class T> class OpMinus
{
 public:
   typedef  T RetType;
   static inline T apply(const T a) {return -a;}
};

//Marker for all trig functions since they return the same units.
template <class T> class OpTrig {typedef T RetType;};

//
//  Trig operators.
//
template <class T> class Opsin
{
 public:
   typedef typename UnaryRetType<T,OpTrig<T> >::RetType RetType;
   static inline RetType apply(const T a) {return sin(a);}
};

template <class T> class Opcos
{
 public:
   typedef typename UnaryRetType<T,OpTrig<T> >::RetType RetType;
   static inline RetType apply(const T a) {return cos(a);}
};

template <class T> class Optan
{
 public:
   typedef typename UnaryRetType<T,OpTrig<T> >::RetType RetType;
   static inline RetType apply(const T a) {return tan(a);}
};

template <class T> class Opasin
{
 public:
   typedef typename UnaryRetType<T,OpTrig<T> >::RetType RetType;
   static inline RetType apply(const T a) {return asin(a);}
};

template <class T> class Opacos
{
 public:
   typedef typename UnaryRetType<T,OpTrig<T> >::RetType RetType;
   static inline RetType apply(const T a) {return acos(a);}
};

template <class T> class Opatan
{
 public:
   typedef typename UnaryRetType<T,OpTrig<T> >::RetType RetType;
   static inline RetType apply(const T a) {return atan(a);}
};

template <class T> class Opsinh
{
 public:
   typedef typename UnaryRetType<T,OpTrig<T> >::RetType RetType;
   static inline RetType apply(const T a) {return sinh(a);}
};

template <class T> class Opcosh
{
 public:
   typedef typename UnaryRetType<T,OpTrig<T> >::RetType RetType;
   static inline RetType apply(const T a) {return cosh(a);}
};

template <class T> class Optanh
{
 public:
   typedef typename UnaryRetType<T,OpTrig<T> >::RetType RetType;
   static inline RetType apply(const T a) {return tanh(a);}
};

//Marker for all exp/log functions since they return the same units.
template <class T> class OpExpLog {typedef T RetType;};

template <class T> class Opexp
{
 public:
   typedef typename UnaryRetType<T,OpExpLog<T> >::RetType RetType;
   static inline RetType apply(const T a) {return exp(a);}
};

template <class T> class Oplog
{
 public:
   typedef typename UnaryRetType<T,OpExpLog<T> >::RetType RetType;
   static inline RetType apply(const T a) {return log(a);}
};

template <class T> class Oplog10
{
 public:
   typedef typename UnaryRetType<T,OpExpLog<T> >::RetType RetType;
   static inline RetType apply(const T a) {return log10(a);}
};

template <class T> class Oppow10
{
 public:
   typedef typename UnaryRetType<T,OpExpLog<T> >::RetType RetType;
   static inline RetType apply(const T a) {return pow10(a);}
};

template <class T> class Oppow2
{
 public:
   typedef typename UnaryRetType<T,Oppow2<T> >::RetType RetType;
   static inline RetType apply(const T a) {return a*a;}
};

template <class T> class Oppow3
{
 public:
   typedef typename UnaryRetType<T,Oppow3<T> >::RetType RetType;
   static inline RetType apply(const T a) {return a*a*a;}
};

template <class T> class Oppow4
{
 public:
   typedef typename UnaryRetType<T,Oppow4<T> >::RetType RetType;
   static inline RetType apply(const T a) {T a2=a*a;return a2*a2;}
};

template <class T> class Opsqrt
{
 public:
   typedef typename UnaryRetType<T,Opsqrt>::RetType RetType;
   static inline RetType apply(const T a) {return sqrt(a);}
};

template <class T> class Opfabs
{
 public:
   typedef  T RetType;
   static inline RetType apply(const T a) {return fabs(a);}
};

template <class T> class Opconj
{
 public:
   typedef  T RetType;
   static inline RetType apply(const T a) {return conj(a);}
};

template <> template <class T> class Opconj<std::complex<T> >
{
 public:
   typedef std::complex<T> RetType;
   static inline std::complex<T> apply(const std::complex<T> a) {return conj(a);}
};

template <class T> class Opabs
{
 public:
   typedef  T RetType;
   static inline RetType apply(const T a) {return abs(a);}
};

template <> template <class T> class Opabs<std::complex<T> >
{
 public:
   typedef  T RetType;
   static inline RetType apply(const std::complex<T> a) {return abs(a);}
};

template <> template <> class Opabs<double>
{
 public:
   typedef  double RetType;
   static inline RetType apply(const double a) {return std::abs(a);}
};

template <class T> class Oparg
{
 public:
   typedef  T RetType;
   static inline RetType apply(const T a) {return arg(a);}
};

template<> template <class T> class Oparg<std::complex<T> >
{
 public:
   typedef  T RetType;
   static inline RetType apply(const std::complex<T> a) {return arg(a);}
};

template <class T> class Opnorm
{
 public:
   typedef  T RetType;
   static inline RetType apply(const T a) {return norm(a);}
};

template<> template <class T> class Opnorm<std::complex<T> >
{
 public:
   typedef  T RetType;
   static inline RetType apply(const std::complex<T> a) {return norm(a);}
};

template <class T> class Opreal
{
 public:
   typedef  T RetType;
   static inline RetType apply(const T a) {return real(a);}
};

template<> template <class T> class Opreal<std::complex<T> >
{
 public:
   typedef  T RetType;
   static inline RetType apply(const std::complex<T> a) {return real(a);}
};

template <class T> class Opimag
{
 public:
   typedef  T RetType;
   static inline RetType apply(const T a) {return imag(a);}
};

template<> template <class T> class Opimag<std::complex<T> >
{
 public:
   typedef  T RetType;
   static inline RetType apply(const std::complex<T> a) {return imag(a);}
};

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

#include "oml/shape.h"
#include "oml/indext.h"
#include "oml/matlimit.h"
//---------------------------------------------------------------
//
//  Temporary unary operation holder.
//
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
  T         operator()(subsc_t n) const {return Op::apply(itsA(n));}
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
  T         operator()(subsc_t i,subsc_t j) const {return Op::apply(itsA(i,j));}
  index_t   size   (                   ) const {return itsA.size();}
  MatLimits GetLimits (                   ) const {return itsA.GetLimits();}
 private:
   A itsA;
};





#endif //_unop_h_
