// File: svbksb.cc

// Modifications for C++ and oml containers Copyright (1994-2003), Jan N. Reimers

#include "oml/vector.h"
#include <cmath>
#include <cassert>
#include <iostream>

#ifndef TYPE
#error "svbksb.cc TYPE not defined"
#endif

template <class T> void SVBackSub(const Matrix<T>& u, const Vector<T>& w,const Matrix<T>& v, const Vector<T>& b, Vector<T>& X)
{
   assert(u.GetColLimits()==w.GetLimits());
   assert(u.GetColLimits()==X.GetLimits());
   assert(u.GetColLimits()==v.GetRowLimits());
   assert(u.GetColLimits()==v.GetColLimits());
   assert(u.GetRowLimits()==b.GetLimits());

   int jj,j,i;
   T    s;
   index_t m=u.GetNumRows();
   index_t n=u.GetNumCols();

   typename Vector<T>::Subscriptor x(X);

   Vector<T> tmp(n);

   for (j=1;j<=n;j++)
   {
     s=0.0;
     if (w(j))
     {
       for (i=1;i<=m;i++) s += u(i,j)*b(i);
       s /= w(j);
     }
     tmp(j)=s;
   }
   for (j=1;j<=n;j++)
   {
     s=0.0;
     for (jj=1;jj<=n;jj++) s += v(j,jj)*tmp(jj);
     x(j)=s;
   }
}

typedef TYPE Type;
template void SVBackSub(const Matrix<Type>&,const Vector<Type>&,
			const Matrix<Type>&,const Vector<Type>&, Vector<Type>&);
