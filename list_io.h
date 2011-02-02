// File: list_io.h  Implementation of an array.
#ifndef _list_io_h_
#define _list_io_h_

// Copyright (1994-2005), Jan N. Reimers

#include "oml/list.h"
#include "oml/array.h"
#include "oml/iterable_io.h"


//-----------------------------------------------------------------------------
//
//  IO
//
template <class T> inline std::ostream& List<T>::Write(std::ostream& os) const
{
  assert(os);
  if (this->Binary()) 
  {
    BinaryWrite(        size(),os);
    BinaryWrite(itsData.size(),os);
  }
  if (this->Ascii())  os << size() << " " << itsData.size() << " ";
//  if (Pretty())  os << "List size=" << size() << ",  Allocated size=" << itsData.size() << std::endl;
  assert(os);
  return ::Write(os,*this);
}

template <class T> inline std::istream& List<T>::Read(std::istream& is)
{
  assert(is);
  index_t next,size;
  if (this->Binary()) 
  {
    BinaryRead(next,is); 
    BinaryRead(size,is);
  }
  if (this->Ascii())
  {
    is >> next >> size;
  }

  assert (next<=size);

  if (this->size()==0)
    SetSize(size);
   else
    assert(itsData.size()==size);

  itsNext=next;

  assert(is);
  return ::Read(is,*this);
}

#endif //_list_io_h_
