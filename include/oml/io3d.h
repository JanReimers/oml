// File: io3d.h  io routines for Vector3D and Matrix3D.
#ifndef _io3d_h_
#define _io3d_h_

// Copyright (1994-2003), Jan N. Reimers

#include "oml/vector3d.h"
#include "oml/matrix3d.h"
#include "oml/imp/stream.h"
#include "oml/imp/binio.h"
#include "oml/imp/mixtypes.h"
#include <iomanip>
#include <cassert>

#include <iosfwd>


//-----------------------------------------------------------------------------
//
//  Vector3D IO.
//
template <class T> std::ostream& operator<<(std::ostream& os,const Vector3D<T>& v)
{
  if (StreamableObject::Binary())
  {
    BinaryWrite(v.x,os);
    BinaryWrite(v.y,os);
    BinaryWrite(v.z,os);
  }
  if (StreamableObject::Ascii()) os << v.x << " " << v.y << " " << v.z << " ";
  if (StreamableObject::Pretty())
  {
    std::streamsize prec=os.precision();
    std::streamsize wid =os.width();
    os << std::setw(0) << "("
       << std::setw(wid) << std::setprecision(prec) << v.x << ","
       << std::setw(wid) << std::setprecision(prec) << v.y << ","
       << std::setw(wid) << std::setprecision(prec) << v.z << ")";
  }
  return os;
}

template <class T> std::istream& operator>>(std::istream& is,Vector3D<T>& v)
{
  if (StreamableObject::Binary())
  {
    BinaryRead(v.x,is);
    BinaryRead(v.y,is);
    BinaryRead(v.z,is);
  }
  if (StreamableObject::Ascii()) is >> v.x >> v.y >> v.z;
  if (StreamableObject::Pretty()) // look for (x,y,z) format.
  {
    char bracket,comma;
    is.get(bracket);
    assert(bracket=='(');
    is >> v.x;
    is.get(comma);
    assert(comma==',');
    is >> v.y;
    is.get(comma);
    assert(comma==',');
    is >> v.z;
    is.get(bracket);
    assert(bracket==')');
  }

  return is;
}

//------------------------------------------------------------------------
//
//  Matrix IO
//
template <class T> std::ostream& operator<<(std::ostream& os,const Matrix3D<T>& a)
{
   if (StreamableObject::Binary())
   {
      BinaryWrite(a.M11,os);
      BinaryWrite(a.M12,os);
      BinaryWrite(a.M13,os);
      BinaryWrite(a.M21,os);
      BinaryWrite(a.M22,os);
      BinaryWrite(a.M23,os);
      BinaryWrite(a.M31,os);
      BinaryWrite(a.M32,os);
      BinaryWrite(a.M33,os);
   }
   if (StreamableObject::Ascii ())
   {
      os << a.M11 << " " << a.M12 << " " << a.M13 << " ";
      os << a.M21 << " " << a.M22 << " " << a.M23 << " ";
      os << a.M31 << " " << a.M32 << " " << a.M33 << " ";
   }
   if (StreamableObject::Pretty())
   {
      std::streamsize prec=os.precision();
      std::streamsize wid =os.width();
      os << std::setw(0);
      os << "[ " << std::setw(wid) << std::setprecision(prec) << a.M11
         << " "  << std::setw(wid) << std::setprecision(prec) << a.M12
         << " "  << std::setw(wid) << std::setprecision(prec) << a.M13
         << " ]" << std::endl;
      os << "[ " << std::setw(wid) << std::setprecision(prec) << a.M21
         << " "  << std::setw(wid) << std::setprecision(prec) << a.M22
         << " "  << std::setw(wid) << std::setprecision(prec) << a.M23
         << " ]" << std::endl;
      os << "[ " << std::setw(wid) << std::setprecision(prec) << a.M31
         << " "  << std::setw(wid) << std::setprecision(prec) << a.M32
         << " "  << std::setw(wid) << std::setprecision(prec) << a.M33
         << " ]" << std::endl;
   }
  return os;
}

template <class T> std::istream& operator>>(std::istream& is,Matrix3D<T>& a)
{
   if (StreamableObject::Binary())
   {
      BinaryRead(a.M11,is);
      BinaryRead(a.M12,is);
      BinaryRead(a.M13,is);
      BinaryRead(a.M21,is);
      BinaryRead(a.M22,is);
      BinaryRead(a.M23,is);
      BinaryRead(a.M31,is);
      BinaryRead(a.M32,is);
      BinaryRead(a.M33,is);
   }
   if (StreamableObject::Ascii ())
   {
      is >> a.M11 >> a.M12 >> a.M13;
      is >> a.M21 >> a.M22 >> a.M23;
      is >> a.M31 >> a.M32 >> a.M33;
   }
   return is;
}


#endif
