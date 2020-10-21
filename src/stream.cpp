// File: StreamableObject.C   Streamable object base class.

// Copyright (1994-2003), Jan N. Reimers

#include "oml/imp/stream.h"
#include <typeinfo>
#include <cassert>
#include <iostream>
#include <string>
#include <cstdlib>

std::string ReadName  (std::istream&);

StreamableObject::Mode StreamableObject::theMode = StreamableObject::ascii;
bool StreamableObject::theCheckArrayTypes=true;

StreamableObject::Mode StreamableObject::SetOutputMode(Mode newMode)
{
  Mode ret=theMode;
  theMode=newMode;
  return ret;
}

bool StreamableObject::CheckArrayTypes(bool check)
{
    bool ret=theCheckArrayTypes;
    theCheckArrayTypes=check;
    return ret;
}

void StreamableObject::WriteHeader(std::ostream& os,c_str type) const
{
  assert(os);
  assert(type);
  if (!Pretty())  os << (int)theMode << " " << type << " ";
  assert(os);
}

StreamableObject::Mode StreamableObject::ReadHeader(std::istream& is,c_str type)
{
  assert(is);
  assert(type);
  Mode current,temp;
  int itemp;
  is >> itemp;
  temp=static_cast<Mode>(itemp);
  assert(temp!=pretty);
  current=SetOutputMode(temp);
  CheckName(is,type);
  assert(is);
  return current;
}

//
//  Check tests if the object type on the stream match this's type.
//  CheckName extracts the name from the stream.
//
void StreamableObject::CheckName(std::istream& is,c_str type)
{
  assert(is);
  assert(type);
  std::string Name=ReadName(is);
  if (theCheckArrayTypes && Name!=type)
  {
    std::cerr << "Read '" << Name << "' from stream, but actual object type is '" << type << "'" << std::endl;
    exit(-1);
  }
}

//
//  ReadName extracts the name from the stream.
//

std::string ReadName(std::istream& is)
{
  assert(is);
  is >> std::ws;
  std::string Name;
  is >>  Name; //Apparently this also reads in the space after the name!
  if (StreamableObject::Binary()) is.get();
  return Name;
}

//
//  Peek at name reads the name from the stream without
//  extracting it.
//
c_str StreamableObject::PeekAtName(std::istream& is)
{

  assert(is);
  is >> std::ws;
  static std::string Name;
  std::streampos beforeName=is.tellg();
  Mode temp;
  int itemp;
  is >> itemp;            //Get the mode flag out of the way.
  (void)temp; //Avoid unused warning
  temp=static_cast<Mode>(itemp);
  is >>  Name;
  is.seekg(beforeName); //Reset stream before start of name and mode flag.
  assert(is);
  return Name.c_str();
}


void OMLArrayIndexError(index_t i, index_t n)
{
  std::cerr << "Array index " << i
	    << " exceeds the array bounds (0:" << n
	    << ")" << std::endl;
}

void OMLListIndexError(index_t i, index_t n)
{
  std::cerr << "List index " << i
	    << " exceeds the list bounds (0:" << n
	    << ")" << std::endl;
}
