// File: tql2.cc

// Modifications for oml containers Copyright (1994-2003), Jan N. Reimers

#include "oml/vector.h"
#include "oml/matrix.h"
#include <cassert>
#include <cmath>
#include <complex>

#ifndef TYPE
#error "tql2.cc TYPE not defined"
#endif

#ifndef MTYPE
#error "tql2.cc MTYPE not defined"
#endif

//----------------------------------------------------------------------
//
//     finds dsqrt(a**2+b**2) without overflow or destructive underflow
//
double pythag(double a,double b)
{
	double p,r,s,t,u;
	double aa=fabs(a),ab=fabs(b);

	p = max(aa,ab);
	if (p == 0.0) return p;

	r = square(min(aa,ab)/p);
	while(true)
	{
		t = 4.0 + r;
		if (t == 4.0) break;
		s = r/t;
		u = 1.0 + 2.0*s;
		p = u*p;
		r = square(s/u) * r;
  }
	return p;
}

template <class T, class M> void tql2(Vector<T>& d, Vector<T>& e, M& z, int& ierr)
{
  int n=d.size();
  assert(d.size()==e.size());
  assert(d.size()==z.GetNumRows());
  assert(d.size()==z.GetNumCols());

  typename Vector<T>::Subscriptor sd(d);
  typename Vector<T>::Subscriptor se(e);
  typename M::Subscriptor sz(z);

  ierr = 0;
  if (n == 1) return;

  for (int j=2;j<=n;j++) se(j-1) = se(j);

  T f = 0.0;
  T tst1 = 0.0;
  se(n) = 0.0;
  for (int l=1;l<=n;l++)
  {
    int j = 0;
    T h = fabs(sd(l)) + fabs(se(l));
    if (tst1 < h) tst1 = h;
    //
    //  look for small sub-diagonal element
    //
    int m;
    for(m=l;m<=n;m++) if (tst1 + fabs(se(m)) == tst1) break;

    if (m!= l)
    {
      do
      {
	if (j== 30)
	{
	  ierr = l;
	  return;
	}
	j=j+1;
	//
	// Form shift.
	//
	int l1 = l + 1;
	int l2 = l1 + 1;
	T g = sd(l);
	T p = (sd(l1) - g) / (2.0 * se(l));
	T r = pythag(p,1.0);
	sd(l) = se(l) / (p + sign(r,p));
	sd(l1) = se(l) * (p + sign(r,p));
	T dl1 = sd(l1);
	h = g - sd(l);
	//		if (l2 <= n)
	for(int k=l2;k<=n;k++) sd(k) = sd(k) - h;
	f = f + h;
	//
	//     .......... ql transformation ..........
	//
	p = sd(m);
	T c = 1.0;
	T c2 = c, c3=c;
	T el1 = se(l1);
	T s = 0.0, s2=0.0;
	//				int mml = m - l;
	for (int ii = 1;ii<=m-l;ii++)
	{

	  c3 = c2;
	  c2 = c;
	  s2 = s;
	  int i = m - ii;
	  g = c * se(i);
	  h = c * p;
	  r = pythag(p,se(i));
	  se(i+1) = s * r;
	  s = se(i) / r;
	  c = p / r;
	  p = c * sd(i) - s * g;
	  sd(i+1) = h + s * (c * g + s * sd(i));
	  for (int k = 1;k<=n;k++) // form vector
	  {
	    h = std::real(sz(k,i+1));
	    sz(k,i+1) = s * sz(k,i) + c * h;
	    sz(k,i) = c * sz(k,i) - s * h;
	  }
	}
	p = -s * s2 * c3 * el1 * se(l) / dl1;
	se(l) = s * p;
	sd(l) = c * p;
      } while (tst1 + fabs(se(l)) > tst1);
    }
    sd(l) = sd(l) + f;
  }
  //
  //  Order eigenvalues and eigenvectors.
  //
  for(int ii = 2;ii<=n;ii++)
  {
    int i = ii - 1;
    int k = i;
    T p = sd(i);

    for(int j = ii;j<=n;j++)
      if (sd(j) < p) {k = j;p = sd(j);}

    if (k!=i)
    {
      sd(k) = sd(i);
      sd(i) = p;
      z.SwapColumns(i,k);
    }
  }
  return;
}

typedef TYPE Type;
typedef MTYPE<std::complex<Type> > MType;

template void tql2<Type,MType>(Vector<Type>& d, Vector<Type>& e, MType& z, int& ierr);

