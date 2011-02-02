// File: binio.h  Inline binario io routines
#ifndef _binio_h_
#define _binio_h_

// Copyright (1994-2005), Jan N. Reimers

#include <iostream>

  template <class T> inline void BinaryWrite(const T& t,std::ostream& os) {os << t;}
  template <class T> inline void BinaryRead (      T& t,std::istream& is) {is >> t;}

#define INTRINSIC_IO(Type) \
  template <> inline void BinaryWrite( Type const &t,std::ostream& os)  {os.write((const char*)&t,sizeof( Type));} \
  template <> inline void BinaryRead ( Type       &t,std::istream& is)  {is.read ((      char*)&t,sizeof( Type));} \

  INTRINSIC_IO(int)
  INTRINSIC_IO(unsigned int)
  INTRINSIC_IO(long int)
  INTRINSIC_IO(unsigned long int)
  INTRINSIC_IO(float)
  INTRINSIC_IO(double)
  INTRINSIC_IO(long double)
  INTRINSIC_IO(short)
  INTRINSIC_IO(char)
  INTRINSIC_IO(bool)
  INTRINSIC_IO(void*)

inline std::istream& operator>>(std::istream& is,void* t) 
{
	long l;
	is >> l >> std::ws ;
	t=reinterpret_cast<void*>(l);
	return is;
} 

#include <complex>

template <typename T> inline void BinaryWrite(const std::complex<T>& t,std::ostream& os) 
{
	BinaryWrite(t.real(),os);
	BinaryWrite(t.imag(),os);
}

template <typename T> inline void BinaryRead(std::complex<T>& t,std::istream& is) 
{
	double re,im;
	BinaryRead(re,is);
	BinaryRead(im,is);
	t=std::complex<T>(re,im);
}

#endif //_binio_h_
