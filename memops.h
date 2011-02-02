// File: memops.h  Template function to implement member op+= type functions.
#ifndef _memops_h_
#define _memops_h_

// Copyright (1994-2005), Jan N. Reimers

#include "oml/indexable.h"
#include <cassert>

#define OP(NAME,OP) \
template <class T, class Derived,Store M,Data D,Shape S> inline \
Derived& Array##NAME (Indexable<T,Derived,M,D,S>& a,const T& scalar)\
{\
  int n=a.size();\
  typename Derived::ArraySubscriptor s(a);\
  for (int i=0;i<n;i++) s[i] OP##=scalar;\
  return static_cast<Derived&>(a);\
}\
template <class T, class Derived,Store M,Data D> inline \
Derived& Vector##NAME (Indexable<T,Derived,M,D,VectorShape>& a,const T& scalar)\
{\
  int hi=a.GetLimits().High;\
  typename Derived::Subscriptor s(a);\
  for (int i=a.GetLimits().Low;i<=hi;i++) s(i) OP##=scalar;\
  return static_cast<Derived&>(a);\
}\
template <class T, class Derived,Store M,Data D> inline \
Derived& Matrix##NAME (Indexable<T,Derived,M,D,MatrixShape>& a,const T& scalar)\
{\
  subsc_t rh=a.GetLimits().Row.High,ch=a.GetLimits().Col.High;\
  typename Derived::Subscriptor s(a);\
  for (int i=a.GetLimits().Row.Low;i<=rh;i++)\
    for (int j=a.GetLimits().Col.Low;j<=ch;j++)\
      s(i,j) OP##=scalar;\
  return static_cast<Derived&>(a);\
}\
template <class T, class Derived,Data D> inline \
Derived& Matrix##NAME (Indexable<T,Derived,Symmetric,D,MatrixShape>& a,const T& scalar)\
{\
  subsc_t rh=a.GetLimits().Row.High,ch=a.GetLimits().Col.High;\
  typename Derived::Subscriptor s(a);\
  for (int i=a.GetLimits().Row.Low;i<=rh;i++)\
    for (int j=i;j<=ch;j++)\
      s(i,j) OP##=scalar;\
  return static_cast<Derived&>(a);\
}\
template <class T,class Derived,Store M,Data D,Shape S,class B> inline \
Derived& Array##NAME (Indexable<T,Derived,M,D,S>& a, const Indexable<T,B,M,D,S>& b)\
{\
	int n=a.size();\
	assert(n==b.size());\
	typename Derived::ArraySubscriptor s(a);\
	for (int i=0;i<n;i++) s[i] OP##= b[i];\
	return static_cast<Derived&>(a);\
}\
template <class T,class Derived,Store M,Data D,class B,Data DB> inline \
Derived& Vector##NAME (Indexable<T,Derived,M,D,VectorShape>& a,\
                  const Indexable<T,B,M,DB,VectorShape>& b)\
{\
  assert(a.GetLimits()==b.GetLimits());\
  int hi=a.GetLimits().High;\
	typename Derived::Subscriptor s(a);\
  for (int i=a.GetLimits().Low;i<=hi;i++) s(i) OP##=b(i);\
	return static_cast<Derived&>(a);\
}\
template <class T,class Derived,Store M,Data D,class B,Data DB> inline \
Derived& Matrix##NAME (Indexable<T,Derived,M,D,MatrixShape>& a,\
                  const Indexable<T,B,M,DB,MatrixShape>& b)\
{\
  assert(a.GetLimits()==b.GetLimits());\
  subsc_t rh=a.GetLimits().Row.High,ch=a.GetLimits().Col.High;\
  typename Derived::Subscriptor s(a);\
  for (int i=a.GetLimits().Row.Low;i<=rh;i++)\
    for (int j=a.GetLimits().Col.Low;j<=ch;j++)\
      s(i,j) OP##=b(i,j);\
  return static_cast<Derived&>(a);\
}\
template <class T,class Derived,Data D,class B,Data DB> inline \
Derived& Matrix##NAME (Indexable<T,Derived,Symmetric,D,MatrixShape>& a,\
                  const Indexable<T,B,Symmetric,DB,MatrixShape>& b)\
{\
  assert(a.GetLimits()==b.GetLimits());\
  subsc_t rh=a.GetLimits().Row.High,ch=a.GetLimits().Col.High;\
  typename Derived::Subscriptor s(a);\
  for (int i=a.GetLimits().Row.Low;i<=rh;i++)\
    for (int j=i;j<=ch;j++)\
      s(i,j) OP##=b(i,j);\
  return static_cast<Derived&>(a);\
}\

OP(Add,+)
OP(Sub,-)
OP(Mul,*)
OP(Div,/)

#undef OP

#ifdef WARN_DEEP_COPY
#include <iostream>
#endif

template <class T, class Derived,Store M,Data D,Shape S,class B> inline
void ArrayAssign(Indexable<T,Derived,M,D,S>& a,const Indexable<T,B,M,D,S>& b)
{
  int n=a.size();
  assert(n==b.size());
  typename Derived::ArraySubscriptor s(a);
  for (int i=0;i<n;i++) s[i]=b[i];
}

template <class T, class Derived,Store M,Data D,class B,Data DB> inline
void VectorAssign(Indexable<T,Derived,M,D,VectorShape>& a,const Indexable<T,B,M,DB,VectorShape>& b)
{
#ifdef WARN_DEEP_COPY
  std::cerr << "Doing abstract VectorAssign n=" << a.size() << std::endl;
#endif
  assert(a.GetLimits()==b.GetLimits());
  int hi=a.GetLimits().High;
  typename Derived::Subscriptor s(a);
  for (int i=a.GetLimits().Low;i<=hi;i++) s(i)=b(i);
}

template <class T, class Derived,Store M,Data D> inline
void VectorAssign(Indexable<T,Derived,M,D,VectorShape>& a,T scalar)
{
#ifdef WARN_DEEP_COPY
  std::cerr << "Doing abstract VectorAssign n=" << a.size() << std::endl;
#endif
  int hi=a.GetLimits().High;
  typename Derived::Subscriptor s(a);
  for (int i=a.GetLimits().Low;i<=hi;i++) s(i)=scalar;
}

template <class T, class Derived, Data D, class B, Store MB, Data DB> inline
void MatrixAssign(Indexable<T,Derived,Full,D,MatrixShape>& a,const Indexable<T,B,MB,DB,MatrixShape>& b)
{
#ifdef WARN_DEEP_COPY
  std::cerr << "Doing abstract MatrixAssign n=" << a.size() << std::endl;
#endif
  assert(a.GetLimits()==b.GetLimits());
  typename Derived::Subscriptor s(a);
  subsc_t rh=a.GetLimits().Row.High,ch=a.GetLimits().Col.High;
  for (int i=a.GetLimits().Row.Low;i<=rh;i++) 
    for (int j=a.GetLimits().Col.Low;j<=ch;j++)
      s(i,j)=b(i,j);
}

template <class T, class Derived, Data D, class B, Store MB, Data DB> inline
void MatrixAssign(Indexable<T,Derived,Symmetric,D,MatrixShape>& a,const Indexable<T,B,MB,DB,MatrixShape>& b)
{
#ifdef WARN_DEEP_COPY
  std::cerr << "Doing abstract MatrixAssign n=" << a.size() << std::endl;
#endif
  assert(a.GetLimits()==b.GetLimits());
  typename Derived::Subscriptor s(a);
  subsc_t rh=a.GetLimits().Row.High,ch=a.GetLimits().Col.High;
  for (int i=a.GetLimits().Row.Low;i<=rh;i++) 
    for (int j=i;j<=ch;j++)
      s(i,j)=b(i,j);
}

#endif //_memops_h_
