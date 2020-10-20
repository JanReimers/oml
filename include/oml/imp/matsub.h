// File: MatrixSub.H  General matrix container functions
#ifndef _MatrixSub_H_
#define _MatrixSub_H_

// Copyright (1994-2003), Jan N. Reimers

#include "oml/imp/index_t.h"
#include "oml/imp/matindex.h"
#include <vector>

template <class T> class Array;
class MatLimits;

template <class T,class A, Store M>
void SetLimits     (Indexable<T,A,M,Real,MatrixShape>& m,const MatLimits& theLimits, bool preserve);
template <class T,class A, Store M>
void ReIndexRows   (Indexable<T,A,M,Real,MatrixShape>& m,const std::vector<index_t>& index);
template <class T,class A, Store M>
void ReIndexColumns(Indexable<T,A,M,Real,MatrixShape>& m,const std::vector<index_t>& index);
template <class T,class A, Store M>
void SwapRows      (Indexable<T,A,M,Real,MatrixShape>& m,index_t i,index_t j);
template <class T,class A, Store M>
void SwapColumns   (Indexable<T,A,M,Real,MatrixShape>& m,index_t i,index_t j);
template <class T,class A, Store M>
void SubMatrix     (Indexable<T,A,M,Real,MatrixShape>& m,const Indexable<T,A,M,Real,MatrixShape>& old);


#endif //_Matrix_H_

