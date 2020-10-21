// File: svdcmp.cc

// Modifications for C++ and oml containers Copyright (1994-2003), Jan N. Reimers

#include "oml/vector.h"
#include "oml/imp/indexsort.h"
#include <cmath>
#include <cassert>
#include <iostream>

#ifndef TYPE
#error "svdcmp.cc TYPE not defined"
#endif

#ifndef MTYPE
#error "svdcmp.cc MTYPE not defined"
#endif


static TYPE at,bt,ct;
#define PYTHAG(a,b) ((at=fabs(a)) > (bt=fabs(b)) ? \
(ct=bt/at,at*sqrt(1.0+ct*ct)) : (bt ? (ct=at/bt,bt*sqrt(1.0+ct*ct)): 0.0))

static TYPE maxarg1,maxarg2;
#define MAX(a,b) (maxarg1=(a),maxarg2=(b),(maxarg1) > (maxarg2) ?\
   (maxarg1) : (maxarg2))
#define SIGN(a,b) ((b) >= 0.0 ? fabs(a) : -fabs(a))

template <class T, class M> void SVDecomp(M& A, Vector<T>& W, M& V)
{
  assert(A.GetColLimits()==W.GetLimits());
  assert(A.GetColLimits()==V.GetRowLimits());
  assert(A.GetColLimits()==V.GetColLimits());

  int flag,l=0,nm=0;
  T c,f,h,s,x,y,z;
  T anorm=0.0,g=0.0,scale=0.0;
  index_t m=A.GetNumRows();
  index_t n=A.GetNumCols();

  if (m < n)
    {
      std::cerr << "SVDCMP: Resizing A with extra zero rows" << std::endl;
      A.SetLimits(n,n,true); //Add the requited rows automatically
      for (int i=m+1;i<=n;i++)
        for (int j=1;j<=n;j++)
            A(i,j)=T(0.0);
    }

  Vector<T> RV1(n);

  typename M::Subscriptor a(A);
  typename M::Subscriptor v(V);
  typename Vector<T>::Subscriptor w(W);
  typename Vector<T>::Subscriptor rv1(RV1);

  for (index_t i=1;i<=n;i++)
    {
      l=i+1;
      rv1(i)=scale*g;
      g=s=scale=0.0;
      if (i <= m)
	{
	  for (index_t k=i;k<=m;k++) scale += fabs(a(k,i));
	  if (scale)
	    {
	      for (index_t k=i;k<=m;k++)
		{
		  a(k,i) /= scale;
		  s += a(k,i)*a(k,i);
		}
	      f=a(i,i);
	      g = -SIGN(sqrt(s),f);
	      h=f*g-s;
	      a(i,i)=f-g;
	      if (i != n)
		{
		  for (index_t j=l;j<=n;j++)
		    {
		      s=0.0;
		      for (index_t k=i;k<=m;k++) s += a(k,i)*a(k,j);
		      f=s/h;
		      for (index_t k=i;k<=m;k++) a(k,j) += f*a(k,i);
		    }
		}
	      for (index_t k=i;k<=m;k++) a(k,i) *= scale;
	    }
	}
      w(i)=scale*g;
      g=s=scale=0.0;
      if (i <= m && i != n)
	{
	  for (index_t k=l;k<=n;k++) scale += fabs(a(i,k));
	  if (scale)
	    {
	      for (index_t k=l;k<=n;k++)
		{
		  a(i,k) /= scale;
		  s += a(i,k)*a(i,k);
		}
	      f=a(i,l);
	      g = -SIGN(sqrt(s),f);
	      h=f*g-s;
	      a(i,l)=f-g;
	      for (index_t k=l;k<=n;k++) rv1(k)=a(i,k)/h;
	      if (i != m) {
		for (index_t j=l;j<=m;j++)
		  {
		    s=0.0;
		    for (index_t k=l;k<=n;k++) s += a(j,k)*a(i,k);
		    for (index_t k=l;k<=n;k++) a(j,k) += s*rv1(k);
		  }
	      }
	      for (index_t k=l;k<=n;k++) a(i,k) *= scale;
	    }
	}
      anorm=MAX(anorm,(fabs(w(i))+fabs(rv1(i))));
    }
  for (index_t i=n;i>=1;i--)
    {
      if (i < n)
	{
	  if (g)
	    {
	      for (index_t j=l;j<=n;j++) v(j,i)=(a(i,j)/a(i,l))/g;
	      for (index_t j=l;j<=n;j++)
		{
		  s=0.0;
		  for (index_t k=l;k<=n;k++) s += a(i,k)*v(k,j);
		  for (index_t k=l;k<=n;k++) v(k,j) += s*v(k,i);
		}
	    }
	  for (index_t j=l;j<=n;j++) v(i,j)=v(j,i)=0.0;
	}
      v(i,i)=1.0;
      g=rv1(i);
      l=i;
    }
  for (index_t i=n;i>=1;i--)
    {
      l=i+1;
      g=w(i);
      if (i < n)
	for (index_t j=l;j<=n;j++) a(i,j)=0.0;
      if (g)
	{
	  g=1.0/g;
	  if (i != n) {
	    for (index_t j=l;j<=n;j++)
	      {
		s=0.0;
		for (index_t k=l;k<=m;k++) s += a(k,i)*a(k,j);
		f=(s/a(i,i))*g;
		for (index_t k=i;k<=m;k++) a(k,j) += f*a(k,i);
	      }
	  }
	  for (index_t j=i;j<=m;j++) a(j,i) *= g;
	}
      else
	{
	  for (index_t j=i;j<=m;j++) a(j,i)=0.0;
	}
      ++a(i,i);
    }
  for (index_t k=n;k>=1;k--)
    {
      for (index_t its=1;its<=30;its++)
	{
	  flag=1;
	  for (l=k;l>=1;l--)
	    {
	      nm=l-1;
	      if (fabs(rv1(l))+anorm == anorm)
		{
		  flag=0;
		  break;
		}
	      if (fabs(w(nm))+anorm == anorm) break;
	    }
	  if (flag)
	    {
	      c=0.0;
	      s=1.0;
	      for (index_t i=l;i<=k;i++)
		{
		  f=s*rv1(i);
		  if (fabs(f)+anorm != anorm)
		    {
		      g=w(i);
		      h=PYTHAG(f,g);
		      w(i)=h;
		      h=1.0/h;
		      c=g*h;
		      s=(-f*h);
		      for (index_t j=1;j<=m;j++)
			{
			  y=a(j,nm);
			  z=a(j,i);
			  a(j,nm)=y*c+z*s;
			  a(j,i)=z*c-y*s;
			}
		    }
		}
	    }
	  z=w(k);
	  if (l == k)
	    {
	      if (z < 0.0)
		{
		  w(k) = -z;
		  for (index_t j=1;j<=n;j++) v(j,k)=(-v(j,k));
		}
	      break;
	    }
	  if (its == 30)
	    {
	      std::cerr << "No convergence in 30 SVDCMP iterations" << std::endl;
	    }
	  x=w(l);
	  nm=k-1;
	  y=w(nm);
	  g=rv1(nm);
	  h=rv1(k);
	  f=((y-z)*(y+z)+(g-h)*(g+h))/(2.0*h*y);
	  g=PYTHAG(f,1.0);
	  f=((x-z)*(x+z)+h*((y/(f+SIGN(g,f)))-h))/x;
	  c=s=1.0;
	  for (index_t j=l;j<=nm;j++)
	    {
	      index_t i=j+1;
	      g=rv1(i);
	      y=w(i);
	      h=s*g;
	      g=c*g;
	      z=PYTHAG(f,h);
	      rv1(j)=z;
	      c=f/z;
	      s=h/z;
	      f=x*c+g*s;
	      g=g*c-x*s;
	      h=y*s;
	      y=y*c;
	      for (index_t jj=1;jj<=n;jj++)
		{
		  x=v(jj,j);
		  z=v(jj,i);
		  v(jj,j)=x*c+z*s;
		  v(jj,i)=z*c-x*s;
		}
	      z=PYTHAG(f,h);
	      w(j)=z;
	      if (z)
		{
		  z=1.0/z;
		  c=f*z;
		  s=h*z;
		}
	      f=(c*g)+(s*y);
	      x=(c*y)-(s*g);
	      for (index_t jj=1;jj<=m;jj++)
		{
		  y=a(jj,j);
		  z=a(jj,i);
		  a(jj,j)=y*c+z*s;
		  a(jj,i)=z*c-y*s;
		}
	    }
	  rv1(l)=0.0;
	  rv1(k)=f;
	  w(k)=x;
	}
    }
    //
    // Sort for descending singular values
    //
    std::vector<index_t> index=MakeDescendingIndex(W); //order of s is preserved.
    W.ReIndex(index);
    A.ReIndexColumns(index);
    V.ReIndexColumns(index);
    //
    //  Reshape to remove temp junk
    //
    if (m < n)
    {
        W.SetLimits(m,true);
        A.SetLimits(m,m,true);
        V.SetLimits(m,n,true);
    }
}

#undef SIGN
#undef MAX
#undef PYTHAG

typedef TYPE Type;
typedef MTYPE<Type> MType;
template void SVDecomp<Type,MType>(MType&,Vector<Type>&,MType&);


