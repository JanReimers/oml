module;
#include <cstddef>
#include <iostream>
#include <cassert>
export module oml.VecIndex;
import oml.VecLimits;
import oml.Indexable;
import oml.ArrIndex;
import oml.Shape;
import oml.Xpr;
import oml.unop;
import oml.StreamableObject;
export
{
// vecindex.h
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
        friend inline bool operator!=(const index_iterator& a, const index_iterator& b) {return a.current!=b.current;}
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

// minmax.h
template <class T> inline       T& Max(      T& a,       T& b)  {return a > b ? a : b;}
template <class T> inline const T& Max(const T& a, const T& b)  {return a > b ? a : b;}
template <class T> inline       T& Min(      T& a,       T& b)  {return a < b ? a : b;}
template <class T> inline const T& Min(const T& a, const T& b)  {return a < b ? a : b;}

// tstram.h
template <class A> class TStreamableObject;

template <class A> std::ostream& operator<<(std::ostream& os, const TStreamableObject<A>& o);
template <class A> std::istream& operator>>(std::istream& os,       TStreamableObject<A>& o);

// Allows template based streaming.  IO methods for the parent type A are called directly.
//  ***No virtual dispatch***!  If you need virtual dispatch see the PMStreamableObject class.
template <class A> class TStreamableObject
: public StreamableObject
{
 public:
  friend std::ostream& operator<< <>(std::ostream& os, const TStreamableObject& o);
  friend std::istream& operator>> <>(std::istream& os,       TStreamableObject& o);
//  friend std::ostream& operator<<(std::ostream& os, const TStreamableObject* o) {return os << *o;}
//  friend std::istream& operator>>(std::istream& is,       TStreamableObject* o) {return is >> *o;}
  static A* Factory(std::istream& is);
};


template <class A> inline std::ostream& operator<<(std::ostream& os, const TStreamableObject<A>& o)
{
  o.WriteHeader(os,typeid(A).name());
  return static_cast<const A&>(o).Write(os);
}

template <class A> inline std::istream& operator>>(std::istream& is,TStreamableObject<A>& o)
{
  StreamableObject::Mode current=o.ReadHeader(is,typeid(A).name());
  static_cast<A&>(o).Read (is);
  o.SetOutputMode(current); //Restore to previous state.
  return is;
}

template <class A> inline A* TStreamableObject<A>::Factory(std::istream& is)
{
  //  CheckName(is,typeid(A).name()); Read operation should do this.
  A* ret=new A;
  is >> ret;
  return ret;
}

}