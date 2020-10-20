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
enum class eOp {Null,Mul};
class OpNull {}; //Marker for operations that don't change dimensions.

template <class TA, class TB> class OpMul; //Forward declare
template <class TA, class TB> class OpDiv; //Forward declare

template <class T1, class T2,class Op=OpNull> struct BinaryRetType;
template <class T1, class T2,eOp=eOp::Null> struct eBinaryRetType;
//
//  Allow case where both types are the same.  Assumes no units.
//
template <class T> struct BinaryRetType<T,T> {typedef T RetType;};
template <class T> struct eBinaryRetType<T,T,eOp::Mul> {typedef T RetType;};
template <class T> struct eBinaryRetType<T,T> {typedef T RetType;};
//
//  complex/scalar mixing is allowed.
//
template <class T,class Op> struct BinaryRetType<             T ,std::complex<T>,Op>
{typedef std::complex<T> RetType;};
template <class T,class Op> struct BinaryRetType<std::complex<T>,             T ,Op>
{typedef std::complex<T> RetType;};
template <class T> struct BinaryRetType<std::complex<T>,std::complex<T>,OpMul<std::complex<T>,std::complex<T> > > {typedef std::complex<T> RetType;};
template <class T> struct BinaryRetType<std::complex<T>,std::complex<T>,OpDiv<std::complex<T>,std::complex<T> > > {typedef std::complex<T> RetType;};

template <class T> struct eBinaryRetType<             T ,std::complex<T>,eOp::Mul> {typedef std::complex<T> RetType;};
template <class T> struct eBinaryRetType<std::complex<T>,             T ,eOp::Mul> {typedef std::complex<T> RetType;};
//template <class T> struct BinaryRetType<std::complex<T>,std::complex<T>,OpMul<std::complex<T>,std::complex<T> > > {typedef std::complex<T> RetType;};
//template <class T> struct BinaryRetType<std::complex<T>,std::complex<T>,OpDiv<std::complex<T>,std::complex<T> > > {typedef std::complex<T> RetType;};
//
//  double int is allowed
//
template <class Op> struct BinaryRetType<double,int   ,Op> {typedef double RetType;};
template <class Op> struct BinaryRetType<int   ,double,Op> {typedef double RetType;};
template <> struct BinaryRetType<double,double,OpMul<double,double> > {typedef double RetType;};
template <> struct BinaryRetType<double,double,OpDiv<double,double> > {typedef double RetType;};
//
//  float int is allowed
//
template <class Op> struct BinaryRetType<float,int  ,Op> {typedef float RetType;};
template <class Op> struct BinaryRetType<int  ,float,Op> {typedef float RetType;};
template <> struct BinaryRetType<float,float,OpMul<float,float> > {typedef float RetType;};
template <> struct BinaryRetType<float,float,OpDiv<float,float> > {typedef float RetType;};
//
//  int int is allowed
//
template <> struct BinaryRetType<int,int,OpMul<int,int> > {typedef int RetType;};
template <> struct BinaryRetType<int,int,OpDiv<int,int> > {typedef int RetType;};


//
//  Primary template for unary operators
//
template <class T,class Op=OpNull> struct UnaryRetType {typedef T RetType;};

template <class T,class Op> struct UnaryRetType<std::complex<T>,Op> {typedef std::complex<T> RetType;};
template <class Op> struct UnaryRetType<double,Op> {typedef double RetType;};
template <class Op> struct UnaryRetType<float ,Op> {typedef float  RetType;};

//
//  Primary template for allowed mxing scalar and array types:
//  Array<Ta>*Ts here we put restictions on Ts only.  The BinaryRetType<Ta,Ts>
//  puts restrictions on Ta.
//

template <class Ts> struct Scalar;

//
//  Now make specializations for each allowed scalar type.
//
#include <complex>
#define ScalarAllowed(S)\
template <> struct Scalar<S>\
{\
  Scalar(S s) : itsVal(s) {};\
  const S& itsVal;\
};\

ScalarAllowed(std::complex<double>)
ScalarAllowed(double)
ScalarAllowed(float)
ScalarAllowed(int)

#undef ScalarAllowed

/*
template <class T> struct Scalar<std::complex<T> >
{
  Scalar(const std::complex<T>& s) : itsVal(s) {};
  const std::complex<T>& itsVal;
};

template <> struct Scalar<std::complex<double> >
{
  Scalar(const std::complex<double>& s) : itsVal(s) {};
  const std::complex<double>& itsVal;
};
*/
#endif // _mixtypes_h_
