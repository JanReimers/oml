// File: minmax.h  Template Min and Max functions.
#ifndef _minmax_h_
#define _minmax_h_

// Copyright (1994-2003), Jan N. Reimers

template <class T> inline       T& Max(      T& a,       T& b)  {return a > b ? a : b;}
template <class T> inline const T& Max(const T& a, const T& b)  {return a > b ? a : b;}
template <class T> inline       T& Min(      T& a,       T& b)  {return a < b ? a : b;}
template <class T> inline const T& Min(const T& a, const T& b)  {return a < b ? a : b;}

#endif
