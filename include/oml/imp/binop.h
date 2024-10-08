// File: binop.h  Glommable Expression Templates.
#ifndef _binop_h_
#define _binop_h_

// Copyright (1994-2003), Jan N. Reimers

#include "oml/imp/mixtypes.h"
#include "oml/imp/shape.h"
#include "oml/imp/index_t.h"
#include "oml/imp/matlimit.h"
#include <cmath>

template <class TR, class TA, class TB, class A, class B, Shape S> class XprBinary
{};


template <class TR, class TA, class TB, class A, class B> class XprBinary<TR,TA,TB,A,B,VectorShape>
{
 public:
   XprBinary(const A& a,const B& b,TR(*f)(const TA&,const TB&)) : itsA(a), itsB(b), itsF(f) {};
  ~XprBinary() {};

  TR        operator[](index_t n) const {return itsF(itsA[n],itsB[n]);}
  TR        operator()(index_t n) const {return itsF(itsA(n),itsB(n));}
  size_t    size      (         ) const {return itsA.size();}
  VecLimits GetLimits (         ) const {return itsA.GetLimits();}
 private:
   A itsA;
   B itsB;
   TR(*itsF)(const TA&,const TB&);
};

template <class TR, class TA, class TB, class A, class B> class XprBinary<TR,TA,TB,A,B,MatrixShape>
{
 public:
   XprBinary(const A& a,const B& b,TR(*f)(const TA&,const TB&)) : itsA(a), itsB(b), itsF(f) {};
  ~XprBinary() {};

  TR         operator[](index_t n          ) const {return itsF(itsA[n],itsB[n]);}
  TR         operator()(index_t i,index_t j) const {return itsF(itsA(i,j),itsB(i,j));}
  size_t     size      (                   ) const {return itsA.size();}
  MatLimits  GetLimits (                   ) const {return itsA.GetLimits();}
 private:
   A itsA;
   B itsB;
   TR(*itsF)(const TA&,const TB&);
};


#endif //_binop_h_
