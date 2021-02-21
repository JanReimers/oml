// File: veclimit.cc   Member functios for VectLimits

// Copyright (1994-2003), Jan N. Reimers

#include "oml/imp/veclimit.h"
#include "oml/imp/stream.h"
#include "oml/imp/binio.h"
#include <iostream>
#include <cassert>




//-----------------------------------------------------------------------------
//
//  IO
//
std::ostream& VecLimits::Write(std::ostream& os) const
{
  assert(os);
  if(StreamableObject::Binary()) {BinaryWrite(Low,os);BinaryWrite(High,os);}
  if(StreamableObject::Ascii ()) os << Low << " " << High << " ";
  if(StreamableObject::Pretty()) os << "(" << Low << ":" << High << ")";
  assert(os);
  return os;
}

std::istream& VecLimits::Read(std::istream& is)
{
  assert(is);
  if(StreamableObject::Binary()) {BinaryRead(Low,is);BinaryRead(High,is);}
  if(StreamableObject::Ascii ()) is >> Low >> High;
  assert(is);
  return is;
}

//-----------------------------------------------------------------------------
//
//  Index checking with descriptive erros.
//
bool VecLimits::Check() const
{
  return High+1 >= Low;
}

bool VecLimits::CheckIndex(index_t i) const
{
  bool ok=true;
  if (i< Low )
  {
    std::cout << "Index " << i << " too low in VecLimits, Low=" << Low << std::endl;
    ok=false;
  }
  if (i>High)
  {
    std::cout << "Index " << i << " too high in VecLimits, High=" << High << std::endl;
    ok=false;
  }
#ifdef DEBUG
  assert(ok);
#endif
  return ok;
}

VecLimits operator*(const VecLimits& a, const VecLimits& b)
{
    assert(a.Low==b.Low);
    int size=a.size()*b.size();
    return VecLimits(a.Low,a.Low+size-1);
}


#undef CHECK
