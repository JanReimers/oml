// File: Ran250.H  Declarations for the R250 random number generator.
#ifndef _Ran250_H_
#define _Ran250_H_

// Copyright (1994-2004), Jan N. Reimers

#include "config.h"
#include <complex>
#include <cassert>
#include <climits>
#ifdef HAVE_STDINT_H
  #include <stdint.h>
#else
  #ifdef HAVE_INTTYPES_H
    #include <inttypes.h>
  #else
    #ifdef  HAVE_SYS_TYPES_H
      #include <sys/types.h>
      typedef  __uint32_t uint32_t;
    #endif
  #endif
#endif

template <class T> union Int2RealUnion;

template <> union Int2RealUnion<float>
{		   	// used to access floats as unsigneds
  float rep;
  ::uint32_t u;
};

template <> union Int2RealUnion<double>
{		   	// used to access doubles as unsigneds
  double rep;
  ::uint32_t u[2];
};

template <class T> class Int2Real;

template <> class Int2Real<float>
{
 public:
  Int2Real();
  float convert(uint32_t l) const;
 private:
  Int2RealUnion<float> itsMantissa;
};

template <> class Int2Real<double>
{
 public:
  Int2Real();
  double convert(uint32_t l1,uint32_t l2) const;
 private:
  Int2RealUnion<double> itsMantissa;
};


inline float Int2Real<float>::convert(uint32_t l) const
{
  Int2RealUnion<float> result;
  result.rep = 1.0;
  result.u   |= (l & itsMantissa.u);
  result.rep -= 1.0;
  assert( result.rep < 1.0 && result.rep >= 0);
  return result.rep;
}

inline double Int2Real<double>::convert(uint32_t l1,uint32_t l2) const
{
  Int2RealUnion<double> result;
  result.rep = 1.0;
  result.u[0] |= (l1 & itsMantissa.u[0]);
  result.u[1] |= (l1 & itsMantissa.u[1]);
  result.rep -= 1.0;
  assert( result.rep < 1.0 && result.rep >= 0);
  return result.rep;
}

class TwoTap
{
  unsigned int   Next,Tap1,Tap2,Mask;
  long*          itsArray;

  Int2Real<float > floatConverter;
  Int2Real<double> doubleConverter;

public:
  TwoTap(unsigned int tap1,unsigned int tap2);
 ~TwoTap();
  const char* Name();

  long GetNext()
  {
    ++Next;
    return itsArray[Next&Mask]=
//      (itsArray[(Next-Tap1)&Mask]%(SHRT_MAX-2))
//       *
//      (itsArray[(Next-Tap2)&Mask]%(SHRT_MAX-2));
      itsArray[(Next-Tap1)&Mask]^itsArray[(Next-Tap2)&Mask];
  }

  float  GetNextFloat () {return  floatConverter.convert(GetNext());}
  double GetNextDouble() {return doubleConverter.convert(GetNext(),GetNext());}


};


class FourTap
{
  unsigned int   Next,Tap1,Tap2,Tap3,Tap4,Mask;
  long*          itsArray;

  Int2Real<float > floatConverter;
  Int2Real<double> doubleConverter;


public:
  FourTap(unsigned int tap1,unsigned int tap2,unsigned int tap3,unsigned int tap4);
 ~FourTap();
  const char* Name();

  long GetNext()
  {
    ++Next;
    return itsArray[Next&Mask]=
      itsArray[(Next-Tap1)&Mask]^
      itsArray[(Next-Tap2)&Mask]^
      itsArray[(Next-Tap3)&Mask]^
      itsArray[(Next-Tap4)&Mask];
  }
  float  GetNextFloat () {return  floatConverter.convert(GetNext());}
  double GetNextDouble() {return doubleConverter.convert(GetNext(),GetNext());}
};

extern FourTap GlobalRandomNumberGenerator;
//extern TwoTap GlobalRandomNumberGenerator;

extern double  INorm;
extern double  LNorm;

template <class T> T      OMLRand();
template <class T> T      OMLRandPos();
template <class T> double OMLRandScale(T max);

template <> inline int    OMLRand<int>   () {return GlobalRandomNumberGenerator.GetNext();}
template <> inline long   OMLRand<long>  () {return GlobalRandomNumberGenerator.GetNext();}
template <> inline float  OMLRand<float> ()
{
  return GlobalRandomNumberGenerator.GetNextFloat();
}
template <> inline double OMLRand<double>()
{
  return GlobalRandomNumberGenerator.GetNextDouble();
}

template <> inline std::complex<double> OMLRand<std::complex<double> >()
{
  return std::complex<double>(OMLRand<double>(),OMLRand<double>());
}

template <> inline int    OMLRandPos<int>   () {return OMLRand<int> ()&0x7fffffff;}
template <> inline long   OMLRandPos<long>  () {return OMLRand<long>()&0x7fffffff;}
template <> inline float  OMLRandPos<float> () {return OMLRand<float>();}
template <> inline double OMLRandPos<double>() {return OMLRand<double>();}
template <> inline std::complex<double> OMLRandPos<std::complex<double> >()
{
  return std::complex<double>(OMLRandPos<double>(),OMLRandPos<double>());
}


template <> inline double OMLRandScale<int>   (int    max) {return INorm*max;}
template <> inline double OMLRandScale<long>  (long   max) {return INorm*max;}
template <> inline double OMLRandScale<float> (float  max) {return max;}
template <> inline double OMLRandScale<double>(double max) {return max;}
template <> inline double OMLRandScale<std::complex<double> >(std::complex<double> max) {return real(max);}

#endif //_Ran250_H_
