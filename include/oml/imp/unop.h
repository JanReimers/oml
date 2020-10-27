// File: unop.h  Glommable Expression Templates.
#ifndef _unop_h_
#define _unop_h_

// Copyright (1994-2005), Jan N. Reimers

#include "oml/imp/shape.h"
#include "oml/imp/index_t.h"
#include "oml/imp/matlimit.h"
#include "oml/imp/mixtypes.h"
#include <cmath>


template <class T> class OpLT
{
 public:
   typedef bool RetType;
   static inline bool apply(const T a, const T b) {return a<b;}
};

template <class T> class OpLE
{
 public:
   typedef bool RetType;
   static inline bool apply(const T a, const T b) {return a<=b;}
};

template <class T> class OpGT
{
 public:
   typedef bool RetType;
   static inline bool apply(const T a, const T b) {return a>b;}
};

template <class T> class OpGE
{
 public:
   typedef bool RetType;
   static inline bool apply(const T a, const T b) {return a>=b;}
};


#endif //_unop_h_
