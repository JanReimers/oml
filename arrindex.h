// File: arrindex.h  Base class for Glommable Expression Templates.
#ifndef _arrindex_h_
#define _arrindex_h_

// Copyright (1994-2005), Jan N. Reimers

#include "oml/indexable.h"
#include "oml/memops.h"

//-------------------------------------------------
//
//  template specialization for Array's.
//
template <class T, class Derived, Store M, Data D> class Indexable<T,Derived,M,D,ArrayShape>
{
 protected: // Users should not be messing around with this class.
  explicit Indexable() {}; 
  Indexable& operator=(const Indexable&) {return *this;}
  Indexable(const Indexable&) {};
  ~Indexable() {};
  
 public:
  T       operator[](index_t n) const {return static_cast<const Derived*>(this)->operator[](n);}
  index_t size      (         ) const {return static_cast<const Derived*>(this)->size();}

  Derived& operator+=(const T& scalar) {return ArrayAdd(*this,scalar);}
  Derived& operator-=(const T& scalar) {return ArraySub(*this,scalar);}
  Derived& operator*=(const T& scalar) {return ArrayMul(*this,scalar);}
  Derived& operator/=(const T& scalar) {return ArrayDiv(*this,scalar);}

  template <class B> Derived& operator+=(const Indexable<T,B,M,D,ArrayShape>& b) {return ArrayAdd(*this,b);}
  template <class B> Derived& operator-=(const Indexable<T,B,M,D,ArrayShape>& b) {return ArraySub(*this,b);}
  template <class B> Derived& operator*=(const Indexable<T,B,M,D,ArrayShape>& b) {return ArrayMul(*this,b);}
  template <class B> Derived& operator/=(const Indexable<T,B,M,D,ArrayShape>& b) {return ArrayDiv(*this,b);}
  
 private:
};


//-------------------------------------------------------------------------
//
//  Some special operators only for arrays.
//
template <class T> class Array;

//
// There is no easy way to get the original derived type out os A which could be
// a hideously complicated expression, so we just create an Array.
//
template <class T, class A, Store M, Data D> inline 
std::ostream& operator<<(std::ostream& os,const Indexable<T,A,M,D,ArrayShape>& a)
{
  return os << Array<T>(a);
}

//--------------------------------------------------------------
//
//  Macros for generating disambiguated template functions that return 
//  expressions.
//
#define DisObOb(Func,Op,SHP)\
template <class TA, class TB, class A, class B, Store M, Data D> inline        \
Xpr<typename Op<TA,TB>::RT,XprBinOp<typename Op<TA,TB>::RT,Ref<TA,Indexable<TA,A,M,D,SHP>,SHP>,  \
Ref<TB,Indexable<TB,B,M,D,SHP>,SHP>,Op<TA,TB>,SHP>,M,D,SHP>                    \
Func (const Indexable<TA,A,M,D,SHP>& a, const Indexable<TB,B,M,D,SHP>& b)      \
{                                                                              \
   typedef Ref<TA,Indexable<TA,A,M,D,SHP>,SHP> RefA;                           \
   typedef Ref<TB,Indexable<TB,B,M,D,SHP>,SHP> RefB;                           \
   typedef XprBinOp<typename Op<TA,TB>::RT,RefA,RefB,Op<TA,TB>,SHP> ExprT;     \
   return Xpr<typename Op<TA,TB>::RT,ExprT,M,D,SHP>(ExprT(RefA(a),RefB(b)));   \
}                                                                              \

DisObOb(operator*,OpMul,ArrayShape)
DisObOb(operator/,OpDiv,ArrayShape)

#endif // _arrindex_h_
