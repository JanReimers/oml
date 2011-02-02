// File: indexable.h  Base class for Glommable Expression Templates.
#ifndef _indexable_h_
#define _indexable_h_

// Copyright (1994-2005), Jan N. Reimers

#include "oml/shape.h"

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
  int n=a.size();
  for (int i=0;i<n;i++) ret+=a[i];
  return ret;
}                                                                                      

template <class A, Store M, Data D, Shape S> inline bool True(const Indexable<bool,A,M,D,S>& a)                                        
{      
  bool ret(true);
  int n=a.size();
  for (int i=0;i<n;i++) ret=ret && a[i];
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
    int n=a.size();
    T ret=n>0 ? a[0] : T(0); // Don't try and read a[0] if there is no data in a!
    for (int i=1;i<n;i++) 
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
  int n=a.size();
  A ret(n);
  for (int i=0;i<n;i++) 
  {
    y0+=a[i];
    ret[i]=y0;
  }
  return ret;
}                                                                                      

template <class T,class A, Store M, Data D, Shape S> inline A Differentiate(const Indexable<T,A,M,D,S>& a)
{      
  int n=a.size();
  A ret(n);
  typename A::ArraySubscriptor s(ret); 
  s[0]=a[0]; //Save integration constant in case caller needs it.
  for (int i=1;i<n;i++) s[i]=a[i]-a[i-1];
  return ret;
}                                                                                      

#include "oml/xpr.h"
#include "oml/binop.h"
#include "oml/unop.h"


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
   typedef XprBinOp<typename Op<Type1,Type2>::RT,RefA,Val<Type2,RefA,S>,Op<Type1,Type2>,S> ExprT;    \
   return Xpr<typename Op<Type1,Type2>::RT,ExprT,M,D,S>(ExprT(RefA(a),Val<Type2,RefA,S>(b,RefA(a))));\
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


//--------------------------------------------------------------
//
//  Macros for generating template functions.
//
#define Ob(Func,Op)                                                                                          \
template <class T, class A, Store M, Data D, Shape S> inline                                                 \
Xpr<typename Op<T>::RetType,XprUnaryOp<typename Op<T>::RetType,Ref<T,Indexable<T,A,M,D,S>,S>,Op<T>,S>,M,D,S> \
Func (const Indexable <T,A,M,D,S>& a)                                                                        \
{                                                                                                            \
  typedef Ref<T,Indexable<T,A,M,D,S>,S> RefA;                                                                \
  typedef XprUnaryOp<typename Op<T>::RetType,RefA,Op<T>,S> ExprT;                                            \
  return Xpr<typename Op<T>::RetType,ExprT,M,D,S>(ExprT(RefA(a)));                                           \
}                                                                                                            \


//----------------------------------------------------------------------------
//
//  Generate lots of expression template functions.
//
ObOb(operator+     ,OpAdd  )
ObOb(operator-     ,OpSub  )
ObOb(DirectMultiply,OpMul  )
ObOb(DirectDivide  ,OpDiv  )
ObOb(operator%     ,OpMod  )
ObOb(pow           ,Oppow  )
ObOb(fmod          ,Opfmod )
ObOb(atan2         ,Opatan2)
ObOb(hypot         ,Ophypot)
ObOb(Equal         ,OpEqual)

ObSc(operator+ ,OpAdd )
ObSc(operator- ,OpSub )
ObSc(operator* ,OpMul )
ObSc(operator/ ,OpDiv )
ObSc(operator% ,OpMod )
ObSc(pow       ,Oppow  )
ObSc(fmod      ,Opfmod )
ObSc(atan2     ,Opatan2)
ObSc(hypot     ,Ophypot)
ObSc(Equal     ,OpEqual)


ObScMix(operator+ ,OpAdd  ,int,double)
ObScMix(operator- ,OpSub  ,int,double)
ObScMix(operator* ,OpMul  ,int,double)
ObScMix(operator/ ,OpDiv  ,int,double)
ObScMix(pow       ,Oppow  ,int,double)
ObScMix(fmod      ,Opfmod ,int,double)
ObScMix(atan2     ,Opatan2,int,double)
ObScMix(hypot     ,Ophypot,int,double)
ObScMix(Equal     ,OpEqual,int,double)

ObScMix(operator+ ,OpAdd  ,std::complex<double>,double)
ObScMix(operator- ,OpSub  ,std::complex<double>,double)
ObScMix(operator* ,OpMul  ,std::complex<double>,double)
ObScMix(operator/ ,OpDiv  ,std::complex<double>,double)
ObScMix(pow       ,Oppow  ,std::complex<double>,double)
ObScMix(fmod      ,Opfmod ,std::complex<double>,double)
ObScMix(atan2     ,Opatan2,std::complex<double>,double)
ObScMix(hypot     ,Ophypot,std::complex<double>,double)
ObScMix(Equal     ,OpEqual,std::complex<double>,double)



  Ob(operator-  , OpMinus  )
  Ob(operator+  , OpPlus  )
  Ob(sin  , Opsin  )
  Ob(cos  , Opcos  )  
  Ob(tan  , Optan  )  
  Ob(asin , Opasin )
  Ob(acos , Opacos )
  Ob(atan , Opatan )
  Ob(sinh , Opsinh )
  Ob(cosh , Opcosh )  
  Ob(tanh , Optanh )  
  Ob(exp  , Opexp  )
  Ob(log  , Oplog  )
  Ob(pow2 , Oppow2 )
  Ob(pow3 , Oppow3 )
  Ob(pow4 , Oppow4 )
  Ob(pow10, Oppow10)
  Ob(log10, Oplog10)
  Ob(sqrt , Opsqrt )
  Ob(fabs , Opfabs )
  Ob(conj , Opconj )
  Ob(abs  , Opabs  )
  Ob(arg  , Oparg  )
  Ob(norm , Opnorm )
  Ob(real , Opreal )
  Ob(imag , Opimag )

template <class T, class A, class B, Store MA, Data DA, Store MB, Data DB, Shape S> inline 
T Dot(const Indexable<T,A,MA,DA,S>& a,const Indexable<T,B,MB,DB,S>& b)
{
  return Sum(DirectMultiply(a,b));
}

template <class T, class A, class B, Store MA, Data DA, Store MB, Data DB, Shape S>
inline std::complex<T> 
Dot(const Indexable<std::complex<T>,A,MA,DA,S>& a,const Indexable<std::complex<T>,B,MB,DB,S>& b)
{
  return Sum(DirectMultiply(a,conj(b)));
}

template <class T, class A, class B, Store MA, Data DA, Store MB, Data DB, Shape S> inline 
typename OpEqual<T,T>::RT operator==(const Indexable<T,A,MA,DA,S>& a, const Indexable<T,B,MB,DB,S>& b)
{
	return True(Equal(a,b));
}

template <class T, class A, class B, Store MA, Data DA, Store MB, Data DB, Shape S> inline 
typename OpEqual<T,T>::RT operator!=(const Indexable<T,A,MA,DA,S>& a, const Indexable<T,B,MB,DB,S>& b)
{
	return !True(Equal(a,b));
}

template <class T, class A, Store M, Data D, Shape S> inline 
typename OpEqual<T,T>::RT operator==(const Indexable<T,A,M,D,S>& a, T b)
{
	return True(Equal(a,b));
}

template <class T, class A, Store M, Data D, Shape S> inline 
typename OpEqual<T,T>::RT operator!=(const Indexable<T,A,M,D,S>& a, T b)
{
	return !True(Equal(a,b));
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
