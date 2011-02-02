// File: unpickle.h  Template functions for unpickling streamable objects.
#ifndef _Unpickle_H_
#define _Unpickle_H_

// Copyright (1994-2003), Jan N. Reimers

#include <fstream>
#include <string>

template <class T> bool UnPickle(T*& pointer, const  char* filep, const char* name)
{
	std::string file(filep);
  bool file_error=true;
  if(file !="") 
  {
    ifstream in(file.c_str(),ios::in | ios::binary);
    if(!in)
    {
      std::cerr << "Can't open " << name << " file :" << file << std::endl;
      file_error=false;
    }
    pointer = T::Factory(in);
    in >> pointer;
  }
  return file_error;
}

//------------------------------------------------------------------
//
//  Non-polymorphic version, so you don't have to define a Factory 
//  function.
//
template <class T> bool UnPickleNP(T*& pointer, const  char* filep, const char* name)
{
  std::string file(filep);
  bool file_error=true;
  if(file !="") 
  {
    ifstream in(file.c_str(),ios::in | ios::binary);
    if(!in)
    {
      std::cerr << "Can't open " << name << " file :" << file << std::endl;
      file_error=false;
    }
    pointer = new T;
    in >> pointer;
  }
  return file_error;
}

#endif // _Unpickle_H_
