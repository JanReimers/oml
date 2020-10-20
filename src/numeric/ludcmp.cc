// File: ludecomp.cc

// Modifications for C++ and oml containers Copyright (1994-2003), Jan N. Reimers

#include "oml/vector.h"
#include <cmath>
#include <cassert>
#include <iostream>

#ifndef TYPE
#error "ludecomp.cc TYPE not defined"
#endif

template <class T> bool LUDecomp(Matrix<T>& A, std::vector<index_t>& ipiv ,T& d)
{
    assert(A.GetRowLimits()==A.GetColLimits());
    assert(A.GetRowLow()==1);
    assert(A.GetColLow()==1);
    T big,dum,sum,temp;
    index_t imax;
    const double TINY=1.0e-20;

    Vector<T> V(A.GetRowLimits());
    index_t n=V.size();

    typename Matrix<T>::Subscriptor a (A);
    typename Vector<T>::Subscriptor vv(V);

    d=1.0;
    for (index_t i=1; i<=n; i++)
    {
        big=0.0;
        for (index_t j=1; j<=n; j++)
            if ((temp=fabs(a(i,j))) > big) big=temp;
        if (big == 0.0)
        {
            std::cerr << "Singular matrix in routine LUDCMP" << std::endl;
            return false;
        }
        vv(i)=1.0/big;
    }
    for (index_t j=1; j<=n; j++)
    {
        for (index_t i=1; i<j; i++)
        {
            sum=a(i,j);
            for (index_t k=1; k<i; k++) sum -= a(i,k)*a(k,j);
            a(i,j)=sum;
        }
        big=-1.0;
        imax=0;
        for (index_t i=j; i<=n; i++)
        {
            sum=a(i,j);
            for (index_t k=1; k<j; k++)
                sum -= a(i,k)*a(k,j);
            a(i,j)=sum;
            if ( (dum=vv(i)*fabs(sum)) >= big)
            {
                big=dum;
                imax=i;
            }
        }
        assert(imax>0);
        if (j != imax)
        {
            for (index_t k=1; k<=n; k++)
            {
                dum=a(imax,k);
                a(imax,k)=a(j,k);
                a(j,k)=dum;
            }
            d = -d;
            vv(imax)=vv(j);
        }
        indx[j-1]=imax;
        if (a(j,j) == 0.0) a(j,j)=TINY;
        if (j != n)
        {
            dum=1.0/(a(j,j));
            for (index_t i=j+1; i<=n; i++) a(i,j) *= dum;
        }
    }
    return true;
}

typedef TYPE Type;
template bool LUDecomp(Matrix<Type>&,Array<index_t>&,Type&);

