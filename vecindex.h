// File: vecindex.h  Base class for Glommable Expression Templates.
#ifndef _vecindex_h_
#define _vecindex_h_

// Copyright (1994-2003), Jan N. Reimers

#include "oml/indexable.h"
#include "oml/memops.h"

//-------------------------------------------------
//
//  template specialization for Vectors's.
//
template <class T, class Derived, Store M, Data D> class Indexable<T,Derived,M,D,VectorShape>
{
 public:

  T operator[](index_t n) const {return static_cast<const Derived*>(this)->operator[](n);}
  T operator()(index_t n) const {return static_cast<const Derived*>(this)->operator()(n);}

  index_t   size  () const {return static_cast<const Derived*>(this)->size();}
  VecLimits GetLimits() const {return static_cast<const Derived*>(this)->GetLimits();}

  Derived& operator+=(const T& scalar) {return ArrayAdd(*this,scalar);}
  Derived& operator-=(const T& scalar) {return ArraySub(*this,scalar);}
  Derived& operator*=(const T& scalar) {return ArrayMul(*this,scalar);}
  Derived& operator/=(const T& scalar) {return ArrayDiv(*this,scalar);}

  template <class B> Derived& operator+=(const Indexable<T,B,M,Real,VectorShape>& b) {return ArrayAdd(*this,b);}
  template <class B> Derived& operator-=(const Indexable<T,B,M,Real,VectorShape>& b) {return ArraySub(*this,b);}
  template <class B> Derived& operator*=(const Indexable<T,B,M,Real,VectorShape>& b) {return ArrayMul(*this,b);}
  template <class B> Derived& operator/=(const Indexable<T,B,M,Real,VectorShape>& b) {return ArrayDiv(*this,b);}

  template <class B> Derived& operator+=(const Indexable<T,B,M,Abstract,VectorShape>& b) {return VectorAdd(*this,b);}
  template <class B> Derived& operator-=(const Indexable<T,B,M,Abstract,VectorShape>& b) {return VectorSub(*this,b);}
  template <class B> Derived& operator*=(const Indexable<T,B,M,Abstract,VectorShape>& b) {return VectorMul(*this,b);}
  template <class B> Derived& operator/=(const Indexable<T,B,M,Abstract,VectorShape>& b) {return VectorDiv(*this,b);}

  class Subscriptor : public Derived::ArraySubscriptor
  {
  public:
    Subscriptor(Indexable& v) : Derived::ArraySubscriptor(v) {};
  };

 protected:
  template <class B> void AssignFrom(const Indexable<T,B,Full,Real    ,VectorShape>& b) { ArrayAssign(*this,b);}
  template <class B> void AssignFrom(const Indexable<T,B,Full,Abstract,VectorShape>& b) {VectorAssign(*this,b);}

  explicit Indexable() {};
  ~Indexable() {};
  Indexable& operator=(const Indexable&) {return *this;}
  Indexable(const Indexable&) {};
};

//--------------------------------------------------------------
//
//  Template specialization for abstract vectors.
//
template <class T, class Derived, Store M> class Indexable<T,Derived,M,Abstract,VectorShape>
{
 public:
  explicit Indexable() {};
  ~Indexable() {};

  T operator()(index_t n) const {return static_cast<const Derived*>(this)->operator()(n);}

  index_t   size  () const {return static_cast<const Derived*>(this)->size();}
  VecLimits GetLimits() const {return static_cast<const Derived*>(this)->GetLimits();}

  Derived& operator+=(T scalar) {return VectorAdd(*this,scalar);}
  Derived& operator-=(T scalar) {return VectorSub(*this,scalar);}
  Derived& operator*=(T scalar) {return VectorMul(*this,scalar);}
  Derived& operator/=(T scalar) {return VectorDiv(*this,scalar);}

  template <class B, Data D> Derived& operator+=(const Indexable<T,B,M,D,VectorShape>& b) {return VectorAdd(*this,b);}
  template <class B, Data D> Derived& operator-=(const Indexable<T,B,M,D,VectorShape>& b) {return VectorSub(*this,b);}
  template <class B, Data D> Derived& operator*=(const Indexable<T,B,M,D,VectorShape>& b) {return VectorMul(*this,b);}
  template <class B, Data D> Derived& operator/=(const Indexable<T,B,M,D,VectorShape>& b) {return VectorDiv(*this,b);}

  class Subscriptor : public Derived::ArraySubscriptor
  {
  public:
    Subscriptor(Indexable& v) : Derived::ArraySubscriptor(v) {};
  };

 protected:
  template <class B,Data DB> void AssignFrom(const Indexable<T,B,M,DB,VectorShape>& b) {VectorAssign(*this,b);}

 private:
  Indexable& operator=(const Indexable&);
  Indexable(const Indexable&);
};



//----------------------------------------------------------------
//
//  Abstract vector specializations for some helper functions.
//

template <class T, class A, Store M> inline
T Sum(const Indexable<T,A,M,Abstract,VectorShape>& a)
{
  T ret(0);
  index_t hi=a.GetLimits().High;
  for (index_t i=a.GetLimits().Low;i<=hi;i++) ret+=a(i);
  return ret;
}

template <class A, Store M> inline
bool True(const Indexable<bool,A,M,Abstract,VectorShape>& a)
{
  bool ret(true);
  index_t hi=a.GetLimits().High;
  for (index_t i=a.GetLimits().Low;i<=hi;i++) ret=ret&&a(i);
  return ret;
}

template <class T, class A, class Op, Store M> class MinMax<T,A,Op,M,Abstract,VectorShape>
{
 public:
  static T apply(const Indexable<T,A,M,Abstract,VectorShape>& a)
  {
    index_t low=a.GetLimits().Low;
    index_t hi =a.GetLimits().High;
    T ret=a(low);
    for (index_t i=low+1;i<=hi;i++)
    {
      T ai=a(i);
      if (Op::apply(ai,ret)) ret=ai;
    }
    return ret;
  }
};

#endif // _vecindex_h_
