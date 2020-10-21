// File: mixtypes.h  Define what types can be mixed by binary operators.
#ifndef _mixtypes_h_
#define _mixtypes_h_

// Copyright (1994-2005), Jan N. Reimers

namespace std {template <class T> class complex;}

//-----------------------------------------------------------
//
//  Specify how different types are allowed to mix in binary operators.
//  Main thing is to define a return type that does not lose any information.
//
template <class T1, class T2> struct ReturnType;
//
//  Allow case where both types are the same.  Assumes no units.
//
template <class T> struct ReturnType<T,T> {typedef T RetType;};
//
//  complex/scalar mixing is allowed.
//
template <class T> struct ReturnType<             T ,std::complex<T>> {typedef std::complex<T> RetType;};
template <class T> struct ReturnType<std::complex<T>,             T > {typedef std::complex<T> RetType;};
//
//  double int is allowed
//
template <> struct ReturnType<double,int   > {typedef double RetType;};
template <> struct ReturnType<int   ,double> {typedef double RetType;};



#endif // _mixtypes_h_
