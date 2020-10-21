// File: cholsky.cc

// Copyright (1994-2003), Jan N. Reimers

#include "oml/matrix.h"
#include "oml/smatrix.h"
#include "oml/numeric.h"
#include <cassert>
#include <cmath>
#include <iostream>
#include <stdlib.h>

#ifndef TYPE
#error "cholsky.cc TYPE not defined"
#endif

//###########################################################################
//
//  Cholsky decomposition of a symmetric-positive definite matrix into
//  upper and lower triangular parts.  A -> U * ~U
//
template <class T> void Cholsky(Matrix<T>& A)
{
    assert(A.GetRowLimits()==A.GetColLimits());
    assert(IsSymmetric(A));
    assert(A.GetRowLow()==1);
    assert(A.GetColLow()==1);

    index_t n=A.GetNumRows();

    typename Matrix<T>::Subscriptor a(A);

    for(index_t j=n; j>=1; j--)
    {
        T temp=0.0;
        for(index_t k=j+1; k<=n; k++) temp+=a(j,k)*a(j,k);
        if (temp>a(j,j))
        {
            std::cerr << "Cholsky(SMatrix<T>& A): Matrix was not positive definite" << std::endl;
            exit(-1);
        }
         a(j,j)=sqrt(a(j,j)-temp);
        for(index_t i=1; i<=j-1; i++)
        {
            temp=0.0;
            for(index_t k=j+1; k<=n; k++) temp+=a(i,k)*a(j,k);
            temp=a(i,j)-temp;
            if (temp!=0) a(i,j)=temp/a(j,j);
            else a(i,j)=0.0;
        }
    }
    for(index_t j=1; j<=n; j++)
        for(index_t i=j+1; i<=n; i++) a(i,j)=0.0;
}

//###########################################################################
//
//  Cholsky decomposition of a symmetric-positive definite matrix into
//  upper and lower triangular parts.  A -> U * ~U
//  Symmetric version, works on upper part, lower part will not be 0.0's.
//
#include "oml/imp/isnan.h"
#define UPPER_ONLY
template <class T> void Cholsky(SMatrix<T>& A)
{
    assert(A.GetRowLow()==1);
    assert(A.GetColLow()==1);

    index_t n=A.GetNumRows();

    typename SMatrix<T>::Subscriptor a(A);

    for(index_t j=n; j>=1; j--)
    {
        T temp=0.0;
        for(index_t k=j+1; k<=n; k++) temp+=a(j,k)*a(j,k);
        if (temp>a(j,j))
        {
            std::cerr << "Cholsky(SMatrix<T>& A): Matrix was not positive definite" << std::endl;
            exit(-1);
        }
        a(j,j)=sqrt(a(j,j)-temp);
        for(index_t i=1; i<=j-1; i++)
        {
            temp=0.0;
            for(index_t k=j+1; k<=n; k++) temp+=a(i,k)*a(j,k);
            temp=a(i,j)-temp;
            if (temp!=0) a(i,j)=temp/a(j,j);
            else a(i,j)=0.0;
        }
    }
}
#undef UPPER_ONLY

typedef TYPE Type;
template void Cholsky(Matrix <Type>&);
template void Cholsky(SMatrix<Type>&);
