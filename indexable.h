// File: indexable.h  Base class for Glommable Expression Templates.
#ifndef _indexable_h_
#define _indexable_h_

// Copyright (1994-2005), Jan N. Reimers

#include "oml/shape.h"
#include "oml/indexable_base.h"
#include <cassert>

namespace std {template <class T> class complex;}

//---------------------------------------------------------------------------
//
//  Primary template for the indexable class.  Use partial specialization
//  for each  container shape.
//
template <class T, class Derived,Store M,Data D,Shape S> class Indexable;

template <class T, class A, Store M, Data D, Shape S> inline T Sum(const Indexable<T,A,M,D,S>& a)
{
  T ret(0);
  for (index_t i: a.array_indices()) ret+=a[i];
  return ret;
}

/*template <class A, Store M, Data D, Shape S> inline bool True(const Indexable<bool,A,M,D,S>& a)
{
  bool ret(true);
  index_t n=a.size();
  for (index_t i=0;i<n;i++) ret=ret && a[i];
  return ret;
}
*/

//------------------------------------------------------------------
//
//  Use classes to workaround lack of support for partial
//  specialization of template functions.
//
template <class T, class A, class Op, Store M, Data D, Shape S> class MinMax
{
public:
    static T apply(const Indexable<T,A,M,D,S>& a)
    {
        T ret=a.size()>0 ? a[0] : T(0); // Don't try and read a[0] if there is no data in a!
        for (index_t i:a.array_indices())
        {
            T ai=a[i];
            if (Op::apply(ai,ret)) ret=ai;
        }
        return ret;
    }
};

//
//  These won't work for expression template, because A won't be a useful return type.
//
template <class T,class A, Store M, Data D, Shape S> inline A Integrate(const Indexable<T,A,M,D,S>& a,T y0=0)
{
  index_t n=a.size();
  A ret(n);
  for (index_t i:a)
  {
    y0+=a[i];
    ret[i]=y0;
  }
  return ret;
}

template <class T,class A, Store M, Data D, Shape S> inline A Differentiate(const Indexable<T,A,M,D,S>& a)
{
  index_t n=a.size();
  A ret(n);
  typename A::ArraySubscriptor s(ret);
  s[0]=a[0]; //Save integration constant in case caller needs it.
  for (index_t i=1;i<n;i++) s[i]=a[i]-a[i-1];
  return ret;
}

#include "oml/xpr.h"
#include "oml/binop.h"
#include "oml/unop.h"

/*
//--------------------------------------------------------------
//
//  Macros for generating template functions that return
//  expressions.
//
#define ObOb(Func,Op)\
template <class TA, class TB, class A, class B, Store M, Data DA, Data DB, Shape S> inline          \
Xpr<typename Op<TA,TB>::RT,XprBinOp<typename Op<TA,TB>::RT,Ref<TA,Indexable<TA,A,M,DA,S>,S>,        \
Ref<TB,Indexable<TB,B,M,DB,S>,S>,Op<TA,TB>,S>,M,BinaryData<DA,DB>::RetData,S>                       \
Func (const Indexable<TA,A,M,DA,S>& a, const Indexable<TB,B,M,DB,S>& b)                             \
{                                                                                                   \
  typedef Ref<TA,Indexable<TA,A,M,DA,S>,S> RefA;                                                    \
  typedef Ref<TB,Indexable<TB,B,M,DB,S>,S> RefB;                                                    \
  typedef XprBinOp<typename Op<TA,TB>::RT,RefA,RefB,Op<TA,TB>,S > ExprT;                            \
  return Xpr<typename Op<TA,TB>::RT,ExprT,M,BinaryData<DA,DB>::RetData,S>(ExprT(RefA(a),RefB(b)));  \
}

#define ObSc(Func,Op)\
template <class T1, class T2, class A, Store M, Data D, Shape S> inline                              \
Xpr<typename Op<T1,T2>::RT,XprBinOp<typename Op<T1,T2>::RT,Ref<T1,Indexable<T1,A,M,D,S>,S>,          \
Val<T2,Ref<T1,Indexable<T1,A,M,D,S>,S>,S>,Op<T1,T2>,S>,M,D,S>                                        \
Func (const Indexable <T1,A,M,D,S>& a, const Scalar<T2>& b)                                          \
{                                                                                                    \
   typedef Ref<T1,Indexable<T1,A,M,D,S>,S> RefA;                                                     \
   typedef XprBinOp<typename Op<T1,T2>::RT,RefA,Val<T2,RefA,S>,Op<T1,T2>,S> ExprT;                   \
   return Xpr<typename Op<T1,T2>::RT,ExprT,M,D,S>(ExprT(RefA(a),Val<T2,RefA,S>(b.itsVal,RefA(a))));  \
}                                                                                                    \
template <class T1, class T2, class B, Store M, Data D, Shape S> inline                              \
Xpr<typename Op<T1,T2>::RT,XprBinOp<typename Op<T1,T2>::RT,Val<T1,Ref<T2,Indexable<T2,B,M,D,S>,S>,S>,\
Ref<T2,Indexable<T2,B,M,D,S>,S>,Op<T1,T2>,S>,M,D,S>                                                  \
Func (const Scalar<T1>& a ,const Indexable <T2,B,M,D,S>& b)                                          \
{                                                                                                    \
   typedef Ref<T2,Indexable<T2,B,M,D,S>,S> RefB;                                                     \
   typedef XprBinOp<typename Op<T1,T2>::RT,Val<T1,RefB,S>,RefB,Op<T1,T2>,S> ExprT;                   \
   return Xpr<typename Op<T1,T2>::RT,ExprT,M,D,S>(ExprT(Val<T1,RefB,S>(a.itsVal,RefB(b)),RefB(b)));  \
}                                                                                                    \
template <class T, class A, Store M, Data D, Shape S> inline                                   \
Xpr<typename Op<T,T>::RT,XprBinOp<typename Op<T,T>::RT,Ref<T,Indexable<T,A,M,D,S>,S>,          \
Val<T,Ref<T,Indexable<T,A,M,D,S>,S>,S>,Op<T,T>,S>,M,D,S>                                       \
Func (const Indexable <T,A,M,D,S>& a, const T& b)                                              \
{                                                                                              \
   typedef Ref<T,Indexable<T,A,M,D,S>,S> RefA;                                                 \
   typedef XprBinOp<typename Op<T,T>::RT,RefA,Val<T,RefA,S>,Op<T,T>,S> ExprT;                  \
   return Xpr<typename Op<T,T>::RT,ExprT,M,D,S>(ExprT(RefA(a),Val<T,RefA,S>(b,RefA(a))));      \
}                                                                                              \
template <class T, class B, Store M, Data D, Shape S> inline                                   \
Xpr<typename Op<T,T>::RT,XprBinOp<typename Op<T,T>::RT,Val<T,Ref<T,Indexable<T,B,M,D,S>,S>,S>, \
Ref<T,Indexable<T,B,M,D,S>,S>,Op<T,T>,S>,M,D,S>                                                \
Func (const T& a ,const Indexable <T,B,M,D,S>& b)                                              \
{                                                                                              \
   typedef Ref<T,Indexable<T,B,M,D,S>,S> RefB;                                                 \
   typedef XprBinOp<typename Op<T,T>::RT,Val<T,RefB,S>,RefB,Op<T,T>,S> ExprT;                  \
   return Xpr<typename Op<T,T>::RT,ExprT,M,D,S>(ExprT(Val<T,RefB,S>(a,RefB(b)),RefB(b)));      \
}

#define ObScMix(Func,Op,Type1,Type2)                                                                 \
template <class A, Store M, Data D, Shape S> inline                                                  \
Xpr<typename Op<Type1,Type2>::RT,XprBinOp<typename Op<Type1,Type2>::RT,                              \
Ref<Type1,Indexable<Type1,A,M,D,S>,S>,                                                               \
Val<Type2,Ref<Type1,Indexable<Type1,A,M,D,S>,S>,S>,Op<Type1,Type2>,S>,M,D,S>                         \
Func (const Indexable <Type1,A,M,D,S>& a, const Type2& b)                                            \
{                                                                                                    \
   typedef Op<Type1,Type2>::RT Type3;\
   typedef Ref<Type1,Indexable<Type1,A,M,D,S>,S> RefA;                                               \
   typedef XprBinOp<Type3,RefA,Val<Type2,RefA,S>,Op<Type1,Type2>,S> ExprT;    \
   return Xpr<Type3,ExprT,M,D,S>(ExprT(RefA(a),Val<Type2,RefA,S>(b,RefA(a))));\
}                                                                                                    \
template <class B, Store M, Data D, Shape S> inline                                                  \
Xpr<typename Op<Type1,Type2>::RT,XprBinOp<typename Op<Type1,Type2>::RT,                              \
Val<Type1,Ref<Type2,Indexable<Type2,B,M,D,S>,S>,S>,                                                  \
Ref<Type2,Indexable<Type2,B,M,D,S>,S>,Op<Type1,Type2>,S>,M,D,S>                                      \
Func (const Type1& a ,const Indexable <Type2,B,M,D,S>& b)                                            \
{                                                                                                    \
   typedef Ref<Type2,Indexable<Type2,B,M,D,S>,S> RefB;                                               \
   typedef XprBinOp<typename Op<Type1,Type2>::RT,Val<Type1,RefB,S>,RefB,Op<Type1,Type2>,S> ExprT;    \
   return Xpr<typename Op<Type1,Type2>::RT,ExprT,M,D,S>(ExprT(Val<Type1,RefB,S>(a,RefB(b)),RefB(b)));\
}                                                                                                    \
template <class A, Store M, Data D, Shape S> inline                                                  \
Xpr<typename Op<Type2,Type1>::RT,XprBinOp<typename Op<Type2,Type1>::RT,                              \
Ref<Type2,Indexable<Type2,A,M,D,S>,S>,                                                               \
Val<Type1,Ref<Type2,Indexable<Type2,A,M,D,S>,S>,S>,Op<Type2,Type1>,S>,M,D,S>                         \
Func (const Indexable <Type2,A,M,D,S>& a, const Type1& b)                                            \
{                                                                                                    \
   typedef Ref<Type2,Indexable<Type2,A,M,D,S>,S> RefA;                                               \
   typedef XprBinOp<typename Op<Type2,Type1>::RT,RefA,Val<Type1,RefA,S>,Op<Type2,Type1>,S> ExprT;    \
   return Xpr<typename Op<Type2,Type1>::RT,ExprT,M,D,S>(ExprT(RefA(a),Val<Type1,RefA,S>(b,RefA(a))));\
}                                                                                                    \
template <class B, Store M, Data D, Shape S> inline                                                  \
Xpr<typename Op<Type2,Type1>::RT,XprBinOp<typename Op<Type2,Type1>::RT,                              \
Val<Type2,Ref<Type1,Indexable<Type1,B,M,D,S>,S>,S>,                                                  \
Ref<Type1,Indexable<Type1,B,M,D,S>,S>,Op<Type2,Type1>,S>,M,D,S>                                      \
Func (const Type2& a ,const Indexable <Type1,B,M,D,S>& b)                                            \
{                                                                                                    \
   typedef Ref<Type1,Indexable<Type1,B,M,D,S>,S> RefB;                                               \
   typedef XprBinOp<typename Op<Type2,Type1>::RT,Val<Type2,RefB,S>,RefB,Op<Type2,Type1>,S> ExprT;    \
   return Xpr<typename Op<Type2,Type1>::RT,ExprT,M,D,S>(ExprT(Val<Type2,RefB,S>(a,RefB(b)),RefB(b)));\
}
*/
//
//  Ob Ob binary
//
template <class TA, class TB, class A, class B, Store M, Data DA, Data DB, Shape S> inline
auto BinaryFunction(const Indexable<TA,A,M,DA,S>& a,
                    const Indexable<TB,B,M,DB,S>& b,
                    typename BinaryRetType<TA,TB>::RetType (*f)(const TA&, const TB&))
{
  typedef typename A::RefT RefA;
  typedef typename B::RefT RefB;
  constexpr Data DR=BinaryData<DA,DB>::RetData ;
  typedef typename BinaryRetType<TA,TB>::RetType TR;
  typedef XprBinaryLambda<TR,TA,TB,RefA,RefB,S> ExprTLambda;
  return Xpr<TR,ExprTLambda,M,DR,S>(ExprTLambda(RefA(a),RefB(b),f));
}

template <class TA, class TB, class A, class B, Store M, Data DA, Data DB, Shape S> inline
auto operator+  (const Indexable<TA,A,M,DA,S>& a, const Indexable<TB,B,M,DB,S>& b)
{
  return BinaryFunction<TA,TB,A,B,M,DA,DB,S>(a,b,[](const TA &xa,const TB& xb) { return xa+xb; });
}

template <class TA, class TB, class A, class B, Store M, Data DA, Data DB, Shape S> inline
auto operator-  (const Indexable<TA,A,M,DA,S>& a, const Indexable<TB,B,M,DB,S>& b)
{
  return BinaryFunction<TA,TB,A,B,M,DA,DB,S>(a,b,[](const TA &xa,const TB& xb) { return xa-xb; });
}

template <class TA, class TB, class A, class B, Store M, Data DA, Data DB, Shape S> inline
auto DirectMultiply (const Indexable<TA,A,M,DA,S>& a, const Indexable<TB,B,M,DB,S>& b)
{
  return BinaryFunction<TA,TB,A,B,M,DA,DB,S>(a,b,[](const TA &xa,const TB& xb) { return xa*xb; });
}
template <class TA, class TB, class A, class B, Store M, Data DA, Data DB, Shape S> inline
auto DirectDivide (const Indexable<TA,A,M,DA,S>& a, const Indexable<TB,B,M,DB,S>& b)
{
  return BinaryFunction<TA,TB,A,B,M,DA,DB,S>(a,b,[](const TA &xa,const TB& xb) { return xa/xb; });
}

//
//  Ob Scalar binary
//
template <class TA, class TB, class A, Store M, Data DA, Shape S> inline
auto BinaryFunction(const Indexable<TA,A,M,DA,S>& a,
                    const TB& b,
                    typename BinaryRetType<TA,TB>::RetType (*f)(const TA&, const TB&))
{
  typedef typename A::RefT RefA;
  typedef  Val<TB,RefA,S> ValB;
  typedef typename BinaryRetType<TA,TB>::RetType TR;
  typedef XprBinaryLambda<TR,TA,TB,RefA,ValB,S> ExprTLambda;
  return Xpr<TR,ExprTLambda,M,DA,S>(ExprTLambda(RefA(a),ValB(b,RefA(a)),f));
}

// Return type is deduced
template <class TA, class TB, class B, Store M, Data DB, Shape S> inline
auto BinaryFunction(const TA& a,
                    const Indexable<TB,B,M,DB,S>& b,
                    typename BinaryRetType<TA,TB>::RetType (*f)(const TA&, const TB&))
{
  typedef typename B::RefT RefB;
  typedef  Val<TA,RefB,S> ValA;
  typedef typename BinaryRetType<TA,TB>::RetType TR;
  typedef XprBinaryLambda<TR,TA,TB,ValA,RefB,S> ExprTLambda;
  return Xpr<TR,ExprTLambda,M,DB,S>(ExprTLambda(ValA(a,RefB(b)),RefB(b),f));
}

// Return type is explicit
template <class TR,class TA, class TB, class A, Store M, Data DA, Shape S> inline
auto BinaryFunctionR(const Indexable<TA,A,M,DA,S>& a,
                    const TB& b,
                    TR (*f)(const TA&, const TB&))
{
  typedef typename A::RefT RefA;
  typedef  Val<TB,RefA,S> ValB;
  typedef XprBinaryLambda<TR,TA,TB,RefA,ValB,S> ExprTLambda;
  return Xpr<TR,ExprTLambda,M,DA,S>(ExprTLambda(RefA(a),ValB(b,RefA(a)),f));
}
template <class TR, class TA, class TB, class B, Store M, Data DB, Shape S> inline
auto BinaryFunctionR(const TA& a,
                    const Indexable<TB,B,M,DB,S>& b,
                    TR (*f)(const TA&, const TB&))
{
  typedef typename B::RefT RefB;
  typedef  Val<TA,RefB,S> ValA;
  typedef XprBinaryLambda<TR,TA,TB,ValA,RefB,S> ExprTLambda;
  return Xpr<TR,ExprTLambda,M,DB,S>(ExprTLambda(ValA(a,RefB(b)),RefB(b),f));
}


#define ObSc1(func,op)\
template <class TA, class A, Store M, Data DA, Shape S> \
inline auto func (const Indexable<TA,A,M,DA,S>& a, const TA& b)\
{  return BinaryFunction<TA,TA,A,M,DA,S>(a,b,[](const TA &xa,const TA& xb) { return op; });}\
template <class TB, class B, Store M, Data DB, Shape S> \
inline auto func  (const TB& a,const Indexable<TB,B,M,DB,S>& b)\
{  return BinaryFunction<TB,TB,B,M,DB,S>(a,b,[](const TB &xa,const TB& xb) { return op; });}



template <class T, class A, class B, class L>
inline bool LogicalII(const Iterable<T,A>& a, const Iterable<T,B>& b,const L& lambda)
{
  assert(a.size()==b.size());
  bool ret(true);
  for (index_t i:a.all())
  {
    ret = ret && lambda(a[i],b[i]);
    if (!ret) break;
  }
  return ret;
}
template <class T, class A, class L>
inline bool Logical(const Iterable<T,A>& a, const T& b,const L& lambda)
{
  bool ret(true);
  for (const T& i:a) ret = ret && lambda(i,b);
  return ret;
}
template <class T, class A, class L>
inline bool Logical(const T & a, const Iterable<T,A>& b,const L& lambda)
{
  bool ret(true);
  for (const T& i:b) ret = ret && lambda(a,i);
  return ret;
}
template <class T, class A, class B,Store MA,Store MB, Data DA, Data DB, class L>
inline bool Logical(const Indexable<T,A,MA,DA,MatrixShape>& a, const Indexable<T,B,MB,DB,MatrixShape>& b,const L& lambda)
{
  assert(a.GetLimits()==b.GetLimits());
  bool ret(true);
    for (index_t i:a.rows())
        for (index_t j:a.cols())
        ret = ret && lambda(a(i,j),b(i,j));
  return ret;
}

#define ObScBool(func,op)\
template <class T, class A,class B,Store M, Data D,Shape S> \
inline bool func(const Indexable<T,A,M,D,S>& a, const Indexable<T,B,M,D,S>& b) \
{return LogicalII(static_cast<const A&>(a),static_cast<const B&>(b),[](const T& xa,const T&xb){return op;});}\
template <class T, class A> inline bool func(const Iterable<T,A>& a, const T& b)\
{return Logical(a,b,[](const T& xa,const T&xb){return op;});}\
template <class T, class A> inline bool func(const T& a, const Iterable<T,A>& b)\
{return Logical(a,b,[](const T& xa,const T&xb){return op;});}\
template <class T, class A,class B,Store MA,Store MB, Data DA, Data DB,Shape S> \
inline bool func(const Indexable<T,A,MA,DA,S>& a, const Indexable<T,B,MB,DB,S>& b) \
{return Logical(a,b,[](const T& xa,const T&xb){return op;});}\

#define ObScMix1(func,op,T1,T2)\
template <class A, Store M, Data DA, Shape S> \
inline auto func (const Indexable<T1,A,M,DA,S>& a, const T2& b) \
{return BinaryFunction<T1,T2,A,M,DA,S>(a,b,[](const T1 &xa,const T2& xb) { return op; });} \
template <class B, Store M, Data DB, Shape S> \
inline auto func  (const T1& a,const Indexable<T2,B,M,DB,S>& b) \
{return BinaryFunction<T1,T2,B,M,DB,S>(a,b,[](const T1 &xa,const T2& xb) { return op; });} \
template <class A, Store M, Data DA, Shape S> \
inline auto func (const Indexable<T2,A,M,DA,S>& a, const T1& b) \
{return BinaryFunction<T2,T1,A,M,DA,S>(a,b,[](const T2 &xa,const T1& xb) { return op; });} \
template <class B, Store M, Data DB, Shape S> \
inline auto func  (const T2& a,const Indexable<T1,B,M,DB,S>& b) \
{return BinaryFunction<T2,T1,B,M,DB,S>(a,b,[](const T2 &xa,const T1& xb) { return op; });}



//----------------------------------------------------------------------------
//
//  Generate lots of expression template functions.
//

ObSc1(operator+ ,xa+xb )
ObSc1(operator- ,xa-xb )
ObSc1(operator* ,xa*xb )
ObSc1(operator/ ,xa/xb )
ObScBool(operator==,xa==xb)
ObScBool(operator!=,xa!=xb)


ObScMix1(operator+ ,xa+xb  ,std::complex<double>,double)
ObScMix1(operator- ,xa-xb  ,std::complex<double>,double)
ObScMix1(operator* ,xa*xb  ,std::complex<double>,double)
ObScMix1(operator/ ,xa/xb  ,std::complex<double>,double)
//ObScMix1(operator+ ,xa+xb  ,std::complex<double>,int)
//ObScMix1(operator- ,xa-xb  ,std::complex<double>,int) //stdcomplex does not support int
//ObScMix1(operator* ,xa*xb  ,std::complex<double>,int)
//ObScMix1(operator/ ,xa/xb  ,std::complex<double>,int)
ObScMix1(operator+ ,xa+xb  ,double,int)
ObScMix1(operator- ,xa-xb  ,double,int)
ObScMix1(operator* ,xa*xb  ,double,int)
ObScMix1(operator/ ,xa/xb  ,double,int)
//ObScMix(Equal     ,OpEqual,std::complex<double>,double)

template <class T, class TR, class A, Store M, Data D, Shape S> inline
auto UnaryFunction(const Indexable <T,A,M,D,S>& a,TR(*f)(const T&))
{
  typedef typename A::RefT RefA;
  typedef XprUnaryLambda<T,TR,RefA,S> ExprTLambda;
  return Xpr<TR,ExprTLambda,M,D,S>(ExprTLambda(RefA(a),f));
}


#define Op(f) \
template <class T, class A, Store M, Data D, Shape S> inline \
auto f(const Indexable <T,A,M,D,S>& a) \
{return UnaryFunction<T,T,A,M,D,S>(a,[](const T &x) { return f(x); });}

#define OpR(f,TR) \
template <class T, class A, Store M, Data D, Shape S> inline \
auto f(const Indexable <T,A,M,D,S>& a) \
{return UnaryFunction<T,TR,A,M,D,S>(a,[](const T &x) { return f(x); });}

template <class T, class A, Store M, Data D, Shape S> inline
auto operator-(const Indexable <T,A,M,D,S>& a)
{return UnaryFunction<T,T,A,M,D,S>(a,[](const T &x) { return -x; });}

template <class T, class A, Store M, Data D, Shape S> inline
auto operator+(const Indexable <T,A,M,D,S>& a)
{return UnaryFunction<T,T,A,M,D,S>(a,[](const T &x) { return x; });}

  Op(sin  )
  Op(cos  )
  Op(tan  )
  Op(asin  )
  Op(acos  )
  Op(atan  )
  Op(sinh  )
  Op(cosh  )
  Op(tanh  )
  Op(exp  )
  Op(log  )
  Op(log10)
  //Op(pow10)
  Op(sqrt)
  Op(conj)
  OpR(real,double)
  OpR(imag,double)
  OpR(norm,double)
  OpR(arg ,double)
  OpR(fabs,double)

  #undef Op

template <class T, class A, class B, Store MA, Data DA, Store MB, Data DB, Shape S> inline
T Dot(const Indexable<T,A,MA,DA,S>& a,const Indexable<T,B,MB,DB,S>& b)
{
  return Sum(DirectMultiply(a,b));
}

/*template <class T, class A, class B, Store MA, Data DA, Store MB, Data DB, Shape S>
inline std::complex<T>
Dot(const Indexable<std::complex<T>,A,MA,DA,S>& a,const Indexable<std::complex<T>,B,MB,DB,S>& b)
{
  return Sum(DirectMultiply(a,conj(b)));
}*/

/*template <class T, class A, class B, Store MA, Data DA, Store MB, Data DB, Shape S> inline
typename OpEqual<T,T>::RT operator==(const Indexable<T,A,MA,DA,S>& a, const Indexable<T,B,MB,DB,S>& b)
{
	return True(Equal(a,b));
}

template <class T, class A, class B, Store MA, Data DA, Store MB, Data DB, Shape S> inline
typename OpEqual<T,T>::RT operator!=(const Indexable<T,A,MA,DA,S>& a, const Indexable<T,B,MB,DB,S>& b)
{
	return !True(Equal(a,b));
}
*/
/*template <class T, class A, Store M, Data D, Shape S> inline
typename OpEqual<T,T>::RT operator==(const Indexable<T,A,M,D,S>& a, T b)
{
	return True(Equal(a,b));
}

template <class T, class A, Store M, Data D, Shape S> inline
typename OpEqual<T,T>::RT operator!=(const Indexable<T,A,M,D,S>& a, T b)
{
	return !True(Equal(a,b));
}
*/
template <class T, class A, Store M, Data D, Shape S> inline T Min(const Indexable<T,A,M,D,S>& a)
{
	return MinMax<T,A,OpLT<T>,M,D,S>::apply(a);
}

template <class T, class A, Store M, Data D, Shape S> inline T Max(const Indexable<T,A,M,D,S>& a)
{
	return MinMax<T,A,OpGT<T>,M,D,S>::apply(a);
}

#endif // _indexable_h_
