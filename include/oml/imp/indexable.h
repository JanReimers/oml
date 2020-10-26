// File: indexable.h  Base class for Glommable Expression Templates.
#ifndef _indexable_h_
#define _indexable_h_

// Copyright (1994-2005), Jan N. Reimers

#include "oml/imp/shape.h"
#include "oml/imp/indexable_base.h"
#include "oml/imp/iterable.h"
#include "oml/imp/xpr.h"
#include "oml/imp/binop.h"
#include "oml/imp/unop.h"
#include <cassert>
#include <complex>

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

template <class T, class A, Data D, Shape S> inline T Sum(const Indexable<T,A,Diagonal,D,S>& a)
{
  T ret(0);
  for (index_t i: a.rows()) ret+=a(i,i);
  return ret;
}
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

//
//  Ob Ob binary
//
template <class TA, class TB, class A, class B, Store MA, Store MB, Data DA, Data DB, Shape S> inline
auto BinaryFunction1(const Indexable<TA,A,MA,DA,S>& a,
                     const Indexable<TB,B,MB,DB,S>& b,
                     typename ReturnType<TA,TB>::RetType (*f)(const TA&, const TB&))
{
  typedef typename A::RefT RefA;
  typedef typename B::RefT RefB;
  constexpr Data  DR=ReturnData <DA,DB>::RetType ;
  constexpr Store MR=ReturnStore<MA,MB>::RetType;
  typedef typename ReturnType<TA,TB>::RetType TR;
  typedef XprBinary<TR,TA,TB,RefA,RefB,S> ExprTLambda;
  return Xpr<TR,ExprTLambda,MR,DR,S>(ExprTLambda(RefA(a),RefB(b),f));
}

template <class TA, class TB, class A, class B, Store MA, Store MB, Data DA, Data DB, Shape S> inline
auto operator+  (const Indexable<TA,A,MA,DA,S>& a, const Indexable<TB,B,MB,DB,S>& b)
{
  return BinaryFunction1<TA,TB,A,B,MA,MB,DA,DB,S>(a,b,[](const TA &xa,const TB& xb) { return xa+xb; });
}

template <class TA, class TB, class A, class B, Store MA, Store MB, Data DA, Data DB, Shape S> inline
auto operator-  (const Indexable<TA,A,MA,DA,S>& a, const Indexable<TB,B,MB,DB,S>& b)
{
  return BinaryFunction1<TA,TB,A,B,MA,MB,DA,DB,S>(a,b,[](const TA &xa,const TB& xb) { return xa-xb; });
}


template <class TA, class TB, class A, class B, Store MA,Store MB, Data DA, Data DB, Shape S> inline
auto DirectMultiply (const Indexable<TA,A,MA,DA,S>& a, const Indexable<TB,B,MB,DB,S>& b)
{
  return BinaryFunction1<TA,TB,A,B,MA,MB,DA,DB,S>(a,b,[](const TA &xa,const TB& xb) { return xa*xb; });
}
template <class TA, class TB, class A, class B, Store MA,Store MB, Data DA, Data DB, Shape S> inline
auto DirectDivide (const Indexable<TA,A,MA,DA,S>& a, const Indexable<TB,B,MB,DB,S>& b)
{
  return BinaryFunction1<TA,TB,A,B,MA,MB,DA,DB,S>(a,b,[](const TA &xa,const TB& xb) { return xa/xb; });
}

//
//  Ob Scalar binary
//
template <class TA, class TB, class A, Store M, Data DA, Shape S> inline
auto BinaryFunction(const Indexable<TA,A,M,DA,S>& a,
                    const TB& b,
                    typename ReturnType<TA,TB>::RetType (*f)(const TA&, const TB&))
{
  typedef typename A::RefT RefA;
  typedef  Val<TB,RefA,S> ValB;
  typedef typename ReturnType<TA,TB>::RetType TR;
  typedef XprBinary<TR,TA,TB,RefA,ValB,S> ExprTLambda;
  return Xpr<TR,ExprTLambda,M,DA,S>(ExprTLambda(RefA(a),ValB(b,RefA(a)),f));
}

// Return type is deduced
template <class TA, class TB, class B, Store M, Data DB, Shape S> inline
auto BinaryFunction(const TA& a,
                    const Indexable<TB,B,M,DB,S>& b,
                    typename ReturnType<TA,TB>::RetType (*f)(const TA&, const TB&))
{
  typedef typename B::RefT RefB;
  typedef  Val<TA,RefB,S> ValA;
  typedef typename ReturnType<TA,TB>::RetType TR;
  typedef XprBinary<TR,TA,TB,ValA,RefB,S> ExprTLambda;
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
  typedef XprBinary<TR,TA,TB,RefA,ValB,S> ExprTLambda;
  return Xpr<TR,ExprTLambda,M,DA,S>(ExprTLambda(RefA(a),ValB(b,RefA(a)),f));
}
template <class TR, class TA, class TB, class B, Store M, Data DB, Shape S> inline
auto BinaryFunctionR(const TA& a,
                    const Indexable<TB,B,M,DB,S>& b,
                    TR (*f)(const TA&, const TB&))
{
  typedef typename B::RefT RefB;
  typedef  Val<TA,RefB,S> ValA;
  typedef XprBinary<TR,TA,TB,ValA,RefB,S> ExprTLambda;
  return Xpr<TR,ExprTLambda,M,DB,S>(ExprTLambda(ValA(a,RefB(b)),RefB(b),f));
}


#define ObSc1(func,op)\
template <class TA, class A, Store M, Data DA, Shape S> \
inline auto func (const Indexable<TA,A,M,DA,S>& a, const TA& b)\
{  return BinaryFunction<TA,TA,A,M,DA,S>(a,b,[](const TA &xa,const TA& xb) { return op; });}\
template <class TB, class B, Store M, Data DB, Shape S> \
inline auto func  (const TB& a,const Indexable<TB,B,M,DB,S>& b)\
{  return BinaryFunction<TB,TB,B,M,DB,S>(a,b,[](const TB &xa,const TB& xb) { return op; });}



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
  typedef XprUnary<T,TR,RefA,S> ExprTLambda;
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

template <class T, class A, Store M, Data D, Shape S> inline T Min(const Indexable<T,A,M,D,S>& a)
{
	return MinMax<T,A,OpLT<T>,M,D,S>::apply(a);
}

template <class T, class A, Store M, Data D, Shape S> inline T Max(const Indexable<T,A,M,D,S>& a)
{
	return MinMax<T,A,OpGT<T>,M,D,S>::apply(a);
}

#endif // _indexable_h_
