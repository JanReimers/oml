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
    for(T i:arr) BinaryRead(i,is);
  else
    for(T& i:arr) is >> i;

  assert(is);
  return is;
}

#endif //_iterable_io_h_

