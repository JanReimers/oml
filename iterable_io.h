// File: iterable_io.h  Implementation of an Iterable object.
#ifndef _iterable_io_h_
#define _iterable_io_h_

// Copyright (1994-2003), Jan N. Reimers

#include "oml/iterable.h"
#include "oml/stream.h"
#include "oml/binio.h"
#include <iostream>
#include <iomanip>
#include <cassert>

//------------------------------------------------------------------------------
//
//  IO stuff.  Just dumps the data to the stream, ... thats it.
//
template <class T, class A> inline std::ostream& Write(std::ostream& os,const Iterable<T,A>& arr)
{
  assert(os);
  typename A::const_iterator b=arr.begin();
  if (StreamableObject::Binary()) for(;b!=arr.end();b++) BinaryWrite(*b,os);
  if (StreamableObject::Ascii ()) for(;b!=arr.end();b++) os << *b << " ";
  if (StreamableObject::Pretty())
  {
    std::streamsize prec=os.precision();
    std::streamsize wid =os.width();
    os << "{ ";
    for(;b!=arr.end();b++) os << std::setw(wid) << std::setprecision(prec) << *b << " ";
    os << "}";
  }
  assert(os);
  return os;
}

template <class T, class A> inline std::istream& Read(std::istream& is,Iterable<T,A>& arr)
{
  assert(is);
  typename A::iterator i=arr.begin();
  if(StreamableObject::Binary())
    for(;i!=arr.end();i++) BinaryRead(*i,is);
  else
    for(;i!=arr.end();i++) is >> *i;

  assert(is);
  return is;
}

#endif //_iterable_io_h_

