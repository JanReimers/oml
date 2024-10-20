// File: vecindex.h  Base class for Glommable Expression Templates.
#ifndef _vecindex_h_
#define _vecindex_h_

// Copyright (1994-2003), Jan N. Reimers

#include "oml/imp/indexable.h"

//
// Template specialization provides index iterators for vector shape
//
template <class Derived> class IndexableBase<Derived,VectorShape>
{
    public:
    //
//  Support range based iteration for rows and columns so client code and do
//     for (index_t i:V)
//          {do something with V(i)
//     for (index_t i:V.all())
//          {do something with V[i]
//
  class index_iterator
  {
    public:
        index_iterator(index_t i) : current{i} {};
        index_iterator operator++(){current++;return (*this);}
        const index_t operator*() const {return current;}
        index_t operator*() {return current;}
        friend bool operator!=(const index_iterator& a, const index_iterator& b) {return a.current!=b.current;}
    private:
        index_t current;
  };

    class iterator_proxy
    {
    public:
        iterator_proxy(const VecLimits& lim) : low(lim.Low), high(lim.High) {};
        iterator_proxy(index_t l, index_t h) : low(l), high(h) {};
        index_iterator begin() const {return low;}
        index_iterator end  () const {return high+1;}
    private:
        index_t low;
        index_t high;
    };

    iterator_proxy indices() const {return iterator_proxy(static_cast<const Derived*>(this)->GetLimits());}
    iterator_proxy indices(index_t i) const 
    {
        return iterator_proxy(i,static_cast<const Derived*>(this)->GetLimits().High);
    }
};

//-------------------------------------------------
//
//  template specialization for Vectors's.
//
template <class T, class Derived, Store M, Data D> class Indexable<T,Derived,M,D,VectorShape>
 : public IndexableBase<Derived,VectorShape>
{
 public:

  T  operator()(index_t n) const {return static_cast<const Derived*>(this)->operator()(n);}
  T& operator()(index_t n)       {return static_cast<      Derived*>(this)->operator()(n);}

  size_t    size  () const {return static_cast<const Derived*>(this)->size();}
  VecLimits GetLimits() const {return static_cast<const Derived*>(this)->GetLimits();}

 protected:
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
 : public IndexableBase<Derived,VectorShape>
{
 public:
  explicit Indexable() {};
  ~Indexable() {};

  T operator()(index_t n) const {return static_cast<const Derived*>(this)->operator()(n);}

  size_t    size  () const {return static_cast<const Derived*>(this)->size();}
  VecLimits GetLimits() const {return static_cast<const Derived*>(this)->GetLimits();}

 private:
  Indexable& operator=(const Indexable&);
  Indexable(const Indexable&);
};

//
//  Create assign functions
//
template <class T, class Derived,Store M,Data D,class B,Data DB> inline
void VectorAssign(Indexable<T,Derived,M,D,VectorShape>& a,const Indexable<T,B,M,DB,VectorShape>& b)
{
#ifdef WARN_DEEP_COPY
  std::cerr << "Doing abstract VectorAssign n=" << a.size() << std::endl;
#endif
  assert(a.GetLimits()==b.GetLimits());
  typename Derived::Subscriptor s(a);
  for (index_t i:a.indices()) s(i)=b(i);
}

template <class T, class Derived,Store M,Data D> inline
void VectorAssign(Indexable<T,Derived,M,D,VectorShape>& a,T scalar)
{
#ifdef WARN_DEEP_COPY
  std::cerr << "Doing abstract VectorAssign n=" << a.size() << std::endl;
#endif
  typename Derived::Subscriptor s(a);
  for (index_t i:a.indices()) s(i)=scalar;
}

#define OP(NAME,OP) \
template <class T, class Derived,Store M,Data D> inline \
Derived& Vector##NAME (Indexable<T,Derived,M,D,VectorShape>& a,const T& scalar)\
{\
  typename Derived::Subscriptor s(a); \
  for (index_t i:a) s(i) OP##=scalar;\
  return static_cast<Derived&>(a);\
}\
template <class T,class Derived,Store M,Data D,class B,Data DB> inline \
Derived& Vector##NAME (Indexable<T,Derived,M,D,VectorShape>& a,\
                  const Indexable<T,B,M,DB,VectorShape>& b)\
{\
    typename Derived::Subscriptor s(a); \
    for (index_t i:a) s(i) OP##=b(i);\
	return static_cast<Derived&>(a);\
}\

OP(Add,+)
OP(Sub,-)
OP(Mul,*)
OP(Div,/)

#undef OP


//----------------------------------------------------------------
//
//  Abstract vector specializations for some helper functions.
//

template <class T, class A, Store M> inline
T Sum(const Indexable<T,A,M,Abstract,VectorShape>& a)
{
  T ret(0);
  for (index_t i:a.indices()) ret+=a(i);
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

template <class T, class A, Store M> inline T Min(const Indexable<T,A,M,Abstract,VectorShape>& a)
{
	return MinMax<T,A,OpLT<T>,M,Abstract,VectorShape>::apply(a);
}

template <class T, class A, Store M> inline T Max(const Indexable<T,A,M,Abstract,VectorShape>& a)
{
	return MinMax<T,A,OpGT<T>,M,Abstract,VectorShape>::apply(a);
}



#endif // _vecindex_h_
