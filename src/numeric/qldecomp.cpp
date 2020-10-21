// File: qldecomp.cc

// Modifications for C++ and oml containers Copyright (1994-2003), Jan N. Reimers

#include "oml/vector.h"
#include <cmath>
#include <cassert>
#include <iostream>

#ifndef TYPE
#error "qldecomp.cc TYPE not defined"
#endif

//##########################################################################
//
//  QL algorithm
//  Calculates Eigenvalues and Eigenstd::vectors of a tridiagonal matrix
//
//  Diagonal     - Input std::vector of diagonal elements
//               - Output std::vector of eigen values
//  OffDiagonal - Input std::vector of Off diagonal elements
//  *this        - Output matrix containing eigen std::vectors as columns.
//
#define SIGN(a,b) ((b)<0 ? -fabs(a) : fabs(a))


template <class T, class M> void QLDecomp(M& A, Vector<T>& Diagonal, Vector<T>& OffDiagonal)
{
  assert(A.GetRowLimits()==A.GetColLimits());
  assert(A.GetRowLow()==1);
  assert(A.GetColLow()==1);

  index_t n=A.GetNumRows();

  T g,r,s,c,p,f,b;

  typename M        ::Subscriptor a (A);
  typename Vector<T>::Subscriptor d (   Diagonal);
  typename Vector<T>::Subscriptor od(OffDiagonal);

  for (index_t i=2;i<=n;i++) od(i-1)=od(i);
  od(n)=0.0;
  for (index_t l=1;l<=n;l++)
  {
    index_t iter=0;
    index_t m;
    do
    {
      for (m=l;m<=n-1;m++)
      {
	T dd=fabs(d(m))+fabs(d(m+1));
	if (fabs(od(m))+dd == dd) break;
      }
      if (m != l)
      {
	if (iter++ == 30)
        {
          std::cerr << "Too many iterations in QL_Decomp";
          assert(false);
        }
	g=(d(l+1)-d(l))/(2.0*od(l));
	r=sqrt((g*g)+1.0);
	g=d(m)-d(l)+od(l)/(g+SIGN(r,g));
	s=c=1.0;
	p=0.0;
	for (index_t i=m-1;i>=l;i--)
	{
	  f=s*od(i);
	  b=c*od(i);
	  if (fabs(f) >= fabs(g))
	  {
	    c=g/f;
	    r=sqrt((c*c)+1.0);
	    od(i+1)=f*r;
	    c *= (s=1.0/r);
	  }
	  else
	  {
	    s=f/g;
	    r=sqrt((s*s)+1.0);
	    od(i+1)=g*r;
	    s *= (c=1.0/r);
	  }
	  g=d(i+1)-p;
	  r=(d(i)-g)*s+2.0*c*b;
	  p=s*r;
	  d(i+1)=g+p;
	  g=c*r-b;
	  for (index_t k=1;k<=n;k++)
	  {
	    f=a(k,i+1);
	    a(k,i+1)=s*a(k,i)+c*f;
	    a(k,i  )=c*a(k,i)-s*f;
	  }
	}
	d(l)=d(l)-p;
	od(l)=g;
	od(m)=0.0;
      }
    } while (m != l);
  }
}

#undef SIGN

typedef TYPE Type;
template void QLDecomp(DMatrix<Type>&,Vector<Type>&,Vector<Type>&);
