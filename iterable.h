// File: Iterable.h  Iterable class for any data type.
#ifndef _Iterable_h_
#define _Iterable_h_

// Copyright (1994-2005), Jan N. Reimers

#include "oml/indext.h"
#include <iosfwd>

// Forward declare some match functions rather than include cmath.
namespace std 
{ 
  long double log(long double);
  long double exp(long double);
}

//-----------------------------------------------------------------------------
//
//  These macros invoke Iterable index bounds checking if DEBUG is on.
//
#if DEBUG
  #include <cassert>
  #define CHECK\
  assert(cur);\
  assert(cur>=start);\
  assert(cur<=end);\
  assert(end>=start);
#else
  #define CHECK
#endif

template <class T, class Derived> class Iterable
{
 protected: //Can only copy and costruct the derived class.
   Iterable() {};
  ~Iterable() {};  
  Iterable(const Iterable&) {};
  Iterable& operator=(const Iterable&) {return *this;}

 public:
  typedef       T*       iterator;
  typedef const T* const_iterator;
  
  const_iterator begin() const {return static_cast<const Derived*>(this)->Get();}
        iterator begin()       {return static_cast<      Derived*>(this)->Get();}
  const_iterator end  () const {return static_cast<const Derived*>(this)->Get()+size();}
        iterator end  ()       {return static_cast<      Derived*>(this)->Get()+size();}

 private:
  index_t size() const {return static_cast<const Derived*>(this)->size();}

};

#undef CHECK


template <class T,class A> inline void Fill(Iterable<T,A>& arr,T value)  
{
  typename A::iterator i=arr.begin();
  for (;i!=arr.end();i++) *i = value;
}

template <class T,class A> inline void FillLinear(Iterable<T,A>& arr,T start, T stop)  
{
  T del = (stop-start)/(double)(static_cast<A&>(arr).size()-1);
  typename A::iterator i=arr.begin();
  for (int n=0;i!=arr.end();i++,n++) *i = start + static_cast<double>(n)*del;
}

template <class A> inline void FillLinear(Iterable<int,A>& arr,int start, int stop)  
{
  int del = (stop-start)/(static_cast<A&>(arr).size()-1);
  typename A::iterator i=arr.begin();
  for (int n=0;i!=arr.end();i++,n++) *i = start + n*del;
}

template <class T,class A> inline void FillPower(Iterable<T,A>& arr,T start, T stop)  
{
  double del=(std::log(stop/start))/(double)(static_cast<A&>(arr).size()-1);
  typename A::iterator i=arr.begin();
  for (int n=0;i!=arr.end();i++,n++) *i=T(start*std::exp(n*del));
}

template <class T1,class T2, template <class T> class A> 
void OML_static_cast(Iterable<T1*,A<T1> >& dest,const Iterable<T2*,A<T2> >& source)  
{
  typename Iterable<T1*,A<T1> >::iterator        i=dest.begin();
  typename Iterable<T2*,A<T2> >::const_iterator  b=source.begin();
  for (;i!=dest.end()&&b!=source.end();i++,b++)  *i = static_cast<T1*>(*b);
}

template <class T,class A> std::ostream& Write  (std::ostream&,const Iterable<T,A>&);
template <class T,class A> std::istream& Read   (std::istream&,      Iterable<T,A>&);



#endif //_Iterable_h_
