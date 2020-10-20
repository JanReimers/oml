// File: array_io.h  inline IO routines for Array's. 
#ifndef _array_io_h_
#define _array_io_h_

// Copyright (1994-2005), Jan N. Reimers

#include "oml/array.h"
#include "oml/iterable_io.h"

//-----------------------------------------------------------------------------
//
//  Write array.
//
template <class T> inline std::ostream& Array<T>::Write(std::ostream& os) const
{
  assert(os);
  if (this->Binary()) BinaryWrite(size(),os);
  if (this->Ascii())  os << size() << " ";
  assert(os);
  return ::Write(os,*this);
}

template <class T> inline std::istream& Array<T>::Read(std::istream& is)
{
  assert(is);
  index_t n;
  if (this->Binary()) BinaryRead(n,is); else is >> n;
  if (size()==0)
    SetSize(n);
   else
    assert(size()==n);

  assert(is);  
  return ::Read(is,*this);
}

#endif //_array_io_h_

