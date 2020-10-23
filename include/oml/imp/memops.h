// File: memops.h  Template function to implement member op+= type functions.
#ifndef _memops_h_
#define _memops_h_

// Copyright (1994-2005), Jan N. Reimers

#include "oml/imp/indexable.h"
#include <cassert>

#define OP(NAME,OP) \
template <class T, class Derived,Store M,Data D,Shape S> inline \
Derived& Array##NAME (Indexable<T,Derived,M,D,S>& a,const T& scalar)\
{\
  typename Derived::ArraySubscriptor s(a); \
  for (index_t i:a.array_indices()) s[i] OP##=scalar;\
  return static_cast<Derived&>(a);\
}\
template <class T, class Derived,Store M,Data D> inline \
Derived& Vector##NAME (Indexable<T,Derived,M,D,VectorShape>& a,const T& scalar)\
{\
  typename Derived::Subscriptor s(a); \
  for (index_t i:a) s(i) OP##=scalar;\
  return static_cast<Derived&>(a);\
}\
template <class T, class Derived,Store M,Data D> inline \
Derived& Matrix##NAME (Indexable<T,Derived,M,D,MatrixShape>& a,const T& scalar)\
{\
  typename Derived::Subscriptor s(a); \
  for (index_t i:a.rows())\
    for (index_t j:a.cols())\
      s(i,j) OP##=scalar;\
  return static_cast<Derived&>(a);\
}\
template <class T, class Derived,Data D> inline \
Derived& Matrix##NAME (Indexable<T,Derived,Symmetric,D,MatrixShape>& a,const T& scalar)\
{\
 typename Derived::Subscriptor s(a); \
 for (index_t i:a.rows())\
    for (index_t j:a.cols())\
      s(i,j) OP##=scalar;\
  return static_cast<Derived&>(a);\
}\
template <class T,class Derived,Store M,Data D,Shape S,class B> inline \
Derived& Array##NAME (Indexable<T,Derived,M,D,S>& a, const Indexable<T,B,M,D,S>& b)\
{\
    typename Derived::ArraySubscriptor s(a); \
	for (index_t i:a.array_indices()) s[i] OP##= b[i];\
	return static_cast<Derived&>(a);\
}\
template <class T,class Derived,Store M,Data D,class B,Data DB> inline \
Derived& Vector##NAME (Indexable<T,Derived,M,D,VectorShape>& a,\
                  const Indexable<T,B,M,DB,VectorShape>& b)\
{\
    typename Derived::Subscriptor s(a); \
    for (index_t i:a) s(i) OP##=b(i);\
	return static_cast<Derived&>(a);\
}\
template <class T,class Derived,Store M,Data D,class B,Data DB> inline \
Derived& Matrix##NAME (Indexable<T,Derived,M,D,MatrixShape>& a,\
                  const Indexable<T,B,M,DB,MatrixShape>& b)\
{\
    typename Derived::Subscriptor s(a); \
    for (index_t i:a.rows())\
        for (index_t j:a.cols())\
            s(i,j) OP##=b(i,j);\
  return static_cast<Derived&>(a);\
}\
template <class T,class Derived,Data D,class B,Data DB> inline \
Derived& Matrix##NAME (Indexable<T,Derived,Symmetric,D,MatrixShape>& a,\
                  const Indexable<T,B,Symmetric,DB,MatrixShape>& b)\
{\
  assert(a.GetLimits()==b.GetLimits());\
  typename Derived::Subscriptor s(a); \
  for (index_t i:a.rows())\
    for (index_t j:a.cols(i))\
      s(i,j) OP##=b(i,j);\
  return static_cast<Derived&>(a);\
}\
/*Diagonal versions*/ \
template <class T,class Derived,Data D,class B,Data DB> inline \
Derived& Matrix##NAME (Indexable<T,Derived,Diagonal,D,MatrixShape>& a,\
                  const Indexable<T,B,Diagonal,DB,MatrixShape>& b)\
{\
    typename Derived::Subscriptor s(a); \
    for (index_t i:a.rows())\
            s(i) OP##=b(i,i);\
  return static_cast<Derived&>(a);\
}\
template <class T, class Derived,Data D> inline \
Derived& Matrix##NAME (Indexable<T,Derived,Diagonal,D,MatrixShape>& a,const T& scalar)\
{\
  typename Derived::Subscriptor s(a); \
  for (index_t i:a.rows())\
      s(i) OP##=scalar;\
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

template <class T, class Derived,Store M,Data DA,Shape S,class B, Data DB> inline
void ArrayAssign(Indexable<T,Derived,M,DA,S>& a,const Indexable<T,B,M,DB,S>& b)
{
  assert(a.size()==b.size());
  for (index_t i:a.array_indices()) a[i]=b[i];
}

template <class T, class Derived,Store M,Data D,class B,Data DB> inline
void VectorAssign(Indexable<T,Derived,M,D,VectorShape>& a,const Indexable<T,B,M,DB,VectorShape>& b)
{
#ifdef WARN_DEEP_COPY
  std::cerr << "Doing abstract VectorAssign n=" << a.size() << std::endl;
#endif
  assert(a.GetLimits()==b.GetLimits());
  typename Derived::Subscriptor s(a);
  for (index_t i:a.indices()) s(i)=b(i);
}

template <class T, class Derived,Store M,Data D> inline
void VectorAssign(Indexable<T,Derived,M,D,VectorShape>& a,T scalar)
{
#ifdef WARN_DEEP_COPY
  std::cerr << "Doing abstract VectorAssign n=" << a.size() << std::endl;
#endif
  typename Derived::Subscriptor s(a);
  for (index_t i:a.indices()) s(i)=scalar;
}

template <class T, class Derived, Data D, class B, Store MB, Data DB> inline
void MatrixAssign(Indexable<T,Derived,Full,D,MatrixShape>& a,const Indexable<T,B,MB,DB,MatrixShape>& b)
{
#ifdef WARN_DEEP_COPY
  std::cerr << "Doing abstract MatrixAssign n=" << a.size() << std::endl;
#endif
  assert(a.GetLimits()==b.GetLimits());
  typename Derived::Subscriptor s(a);
  #pragma omp parallel for collapse(2)
  for (index_t i=a.GetLimits().Row.Low;i<=a.GetLimits().Row.High;i++)
    for (index_t j=a.GetLimits().Col.Low;j<=a.GetLimits().Col.High;j++)
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
  for (index_t i:a.rows())
    for (index_t j:a.rows(i))
      s(i,j)=b(i,j);
}

template <class T, class Derived, Data D, class B, Data DB> inline
void DiagonalAssign(Indexable<T,Derived,Diagonal,D,MatrixShape>& a,const Indexable<T,B,Diagonal,DB,MatrixShape>& b)
{
#ifdef WARN_DEEP_COPY
  std::cerr << "Doing abstract MatrixAssign n=" << a.size() << std::endl;
#endif
  assert(a.GetLimits()==b.GetLimits());
  typename Derived::Subscriptor s(a);
  for (index_t i:a.rows())
      s(i)=b(i,i);
}

#endif //_memops_h_
