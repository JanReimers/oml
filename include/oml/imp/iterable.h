// File: Iterable.h  Iterable class for any data type.
#ifndef _Iterable_h_
#define _Iterable_h_

// Copyright (1994-2005), Jan N. Reimers

#include "oml/imp/index_t.h"
#include "oml/imp/stream.h"
#include "oml/imp/binio.h"
#include <iosfwd>
#include <iomanip>
#include <cmath>
#include <complex>

//-----------------------------------------------------------------------------
//
//  These macros invoke Iterable index bounds checking if DEBUG is on.
//
#if DEBUG
  #include <cassert>
  #define CHECK(i)\
  assert(i>=0);\
  assert(i<size());
#else
  #define CHECK(i)
#endif

template <class T, class Derived> class Iterable
{
protected: //Can only copy and construct the derived class.
    Iterable() {};
    ~Iterable() {};
    Iterable(const Iterable&) {};
    Iterable& operator=(const Iterable&) {return *this;}

public:
    typedef       T*       iterator;
    typedef const T* const_iterator;

    T operator[](index_t i) const {CHECK(i);return begin()[i];}

    const_iterator begin() const {return static_cast<const Derived*>(this)->priv_begin();}
          iterator begin()       {return static_cast<      Derived*>(this)->priv_begin();}
    const_iterator end  () const {return static_cast<const Derived*>(this)->priv_begin()+size();}
          iterator end  ()       {return static_cast<      Derived*>(this)->priv_begin()+size();}

    index_t size() const {return static_cast<const Derived*>(this)->size();}

    //
    //  All for (index_t i:a.indices()) range loops over indices
    //  This is useful when you want to iterate over two containers.
    //
    class index_iterator
    {
    public:
        index_iterator(index_t i) : current{i} {};
        index_iterator operator++(){current++;return (*this);}
        const index_t operator*() const {return current;}
              index_t operator*() {return current;}
        bool operator!=(const index_iterator& b) {return current!=b.current;}
    private:
        index_t current;
    };

    class iterator_proxy
    {
    public:
        iterator_proxy(index_t l, index_t h) : low(l), high(h) {};
        index_iterator begin() const {return low;}
        index_iterator end  () const {return high+1;}
    private:
        index_t low;
        index_t high;
    };

    iterator_proxy indices() const {return iterator_proxy(0,size()-1);}

};

#undef CHECK


template <class T,class A> inline void Fill(Iterable<T,A>& arr,T value)
{
  for (T& a:arr) a = value; //Sill fails for complex
}

template <class T,class A> inline void FillLinear(Iterable<T,A>& arr,T start, T stop)
{
  T del = (stop-start)/(double)(static_cast<A&>(arr).size()-1);
  T val=start;
  for (T& a:arr)
  {
      a = val;
      val+=del;
  }
}

//------------------------------------------------------------------------------
//
//  IO stuff.  Just dumps the data to the stream, ... thats it.
//
template <class T, class A> inline std::ostream& Write(std::ostream& os,const Iterable<T,A>& arr)
{
  assert(os);
  if (StreamableObject::Binary()) for(T b:arr) BinaryWrite(b,os);
  if (StreamableObject::Ascii ()) for(T b:arr) os << b << " ";
  if (StreamableObject::Pretty())
  {
    std::streamsize prec=os.precision();
    std::streamsize wid =os.width();
    os << "{ ";
    for(T b:arr) os << std::setw(wid) << std::setprecision(prec) << b << " ";
    os << "}";
  }
  assert(os);
  return os;
}

template <class T, class A> inline std::istream& Read(std::istream& is,Iterable<T,A>& arr)
{
  assert(is);
  if(StreamableObject::Binary())
    for(T& i:arr) BinaryRead(i,is);
  else
    for(T& i:arr) is >> i;

  assert(is);
  return is;
}


//
//  Logical operators mapped over iterable arrays
//
template <class T, class A, class B, class L>
inline bool LogicalII(const Iterable<T,A>& a, const Iterable<T,B>& b,const L& lambda)
{
  assert(a.size()==b.size());
  bool ret(true);
  for (index_t i:a.indices())
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

//inline constexpr bool isnan(const std::complex<double>& c)
//{
//    return std::isnan(c.real()) || std::isnan(c.imag());
//}
inline constexpr bool isnan(const std::complex<double>& c)
{
    return std::isnan(c.real()) || std::isnan(c.imag());
}
inline bool isinf(const std::complex<double>& c)
{
    return std::isinf(c.real()) || std::isinf(c.imag());
}

template <class T, class D> inline bool isnan(const Iterable<T,D>& arr)
{
    bool ret=false;
    for (const T& a : arr) if (isnan(a)) {ret=true;break;}
    return ret;
}

template <class T, class D> inline bool isinf(const Iterable<T,D>& arr)
{
    bool ret=false;
    for (const T& a : arr) if (isinf(a)) {ret=true;break;}
    return ret;
}




#endif //_Iterable_h_
