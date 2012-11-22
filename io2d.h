// File: io2d.h  io routines for Vector2d and Matrix2d.
#ifndef _io2d_h_
#define _io2d_h_

// Copyright (1994-2003), Jan N. Reimers

#include "oml/vector2d.h"
#include "oml/stream.h"
#include "oml/binio.h"
#include "oml/mixtypes.h"
#include <iomanip>
#include <cassert>

#include <iosfwd>


//-----------------------------------------------------------------------------
//
//  Vector2d IO.
//
template <class T> std::ostream& operator<<(std::ostream& os,const Vector2d<T>& v)
{
  if (StreamableObject::Binary())
  {
    BinaryWrite(v.x,os);
    BinaryWrite(v.y,os);
  }
  if (StreamableObject::Ascii()) os << v.x << " " << v.y << " ";
  if (StreamableObject::Pretty())
  {
    int prec=os.precision();
    int wid =os.width();
    os << std::setw(0) << "("
       << std::setw(wid) << std::setprecision(prec) << v.x << ","
       << std::setw(wid) << std::setprecision(prec) << v.y << ")";
  }
  return os;
}

template <class T> std::istream& operator>>(std::istream& is,Vector2d<T>& v)
{
  if (StreamableObject::Binary())
  {
    BinaryRead(v.x,is);
    BinaryRead(v.y,is);
  }
  if (StreamableObject::Ascii()) is >> v.x >> v.y;
  if (StreamableObject::Pretty()) // look for (x,y) format.
  {
    char bracket,comma;
    is.get(bracket);
    assert(bracket=='(');
    is >> v.x;
    is.get(comma);
    assert(comma==',');
    is >> v.y;
    is.get(bracket);
    assert(bracket==')');
  }

  return is;
}




#endif
