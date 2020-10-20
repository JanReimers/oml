// File: Random.h  Functions for filling an array with random numbers.
#ifndef _Random_h_
#define _Random_h_

// Copyright (1994-2003), Jan N. Reimers

#include "oml/iterable.h"
#include "oml/ran250.h"

/*! \file random.h
  \brief Routines for filling OML containers with random numbers.

  Random number are supplied by a 2 tap xor random number generator. Very fast!
 */
//! Fill with random numbers, range = (-1,1).
template <class T,class A> void FillRandom(Iterable<T,A>& arr)
{
  for (T& i:arr) i=OMLRand<T>();
}

//! Fill with positive random numbers, range = [0,1).
template <class T,class A> void FillRandomPositive(Iterable<T,A>& arr)
{
  for (T& i:arr) i=OMLRandPos<T>();
}

//! Fill with positive random numbers, range = [-max,max).
template <class T,class A> void FillRandom(Iterable<T,A>& arr,T max)
{
  double scale=OMLRandScale<T>(max);
  for (T& i:arr) i = (T)(OMLRand<T>() * scale);
}

//! Fill with positive random numbers, range = [0,max).
template <class T,class A> void FillRandomPositive(Iterable<T,A>& arr,T max)
{
  double scale=OMLRandScale<T>(max);
  for (T& i:arr) i = (T)(OMLRandPos<T>() * scale);
}

#endif //_Random_h_
