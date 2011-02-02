// File: TStreamableObject.H  Template mixin class for streamable objects.
#ifndef _TStreamableObject_H_
#define _TStreamableObject_H_

// Copyright (1994-2005), Jan N. Reimers

#include "oml/stream.h"
#include <typeinfo>

template <class A> class TStreamableObject;

template <class A> std::ostream& operator<<(std::ostream& os, const TStreamableObject<A>& o);
template <class A> std::istream& operator>>(std::istream& os,       TStreamableObject<A>& o);

template <class A> class TStreamableObject
: public StreamableObject
{
 public:
  friend std::ostream& operator<< <>(std::ostream& os, const TStreamableObject& o);
  friend std::istream& operator>> <>(std::istream& os,       TStreamableObject& o);
  friend std::ostream& operator<<(std::ostream& os, const TStreamableObject* o) {return os << *o;}
  friend std::istream& operator>>(std::istream& is,       TStreamableObject* o) {return is >> *o;}
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

#endif //_TStreamabelObject_H_
