// File: matlimit.C   Member functions for MatLimits.

// Copyright (1994-2003), Jan N. Reimers

#include "oml/matlimit.h"
#include "oml/stream.h"
#include <iostream>
#include <cassert>

#ifdef DEBUG
#define CHECK Check()
#else
#define CHECK
#endif


//------------------------------------------------------------------------
//
//  Checking functions.
//
bool MatLimits::CheckIndex(index_t row, index_t col) const
{
  bool ok=true;
  if (!Row.CheckIndex(row))
  {
    std::cout << "Matrix row index " << row << " out of range " << Row << std::endl;
    ok=false;
  }
  if (!Col.CheckIndex(col))
  {
    std::cout << "Matrix column index " << col << " out of range " << Col << std::endl;
    ok=false;
  }
#ifdef DEBUG
  assert(ok);
#endif
  return ok;
}


//------------------------------------------------------------------------
//
// IO
//
std::ostream& MatLimits::Write(std::ostream& os) const
{
  assert(os);
  if(StreamableObject::Binary()) os << Row << Col;
  if(StreamableObject::Ascii ()) os << Row << " " << Col << " ";
  if(StreamableObject::Pretty()) os << Row << "," << Col << " ";
  assert(os);
  return os;
}
std::istream& MatLimits::Read(std::istream& is)
{
  assert(is);
  if(StreamableObject::Binary()) is >> Row >> Col;
  if(StreamableObject::Ascii ()) is >> Row >> Col;
  assert(is);
  return is;
}
