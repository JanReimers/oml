// File: tgl2.cc

// Modifications for oml containers Copyright (1994-2003), Jan N. Reimers

#include "oml/vector.h"
#include "oml/matrix.h"
#include <cassert>
#include <cmath>
#include <complex>

inline double max(double a, double b) {return a>b ? a : b;}
inline double min(double a, double b) {return a<b ? a : b;}
inline double square(double a) {return a*a;}
inline double sign(double a, double b) {return b<0 ? -fabs(a) : fabs(a);}


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

void tql2(Vector<double>& d, Vector<double>& e, Matrix<std::complex<double> >& z, int ierr)
{
  int n=d.size();
  assert(d.size()==e.size());
  assert(d.size()==z.GetNumRows());
  assert(d.size()==z.GetNumCols());

	Vector<double>::Subscriptor sd(d);
	Vector<double>::Subscriptor se(e);
	Matrix<std::complex<double> >::Subscriptor sz(z);

  ierr = 0;
  if (n == 1) return;

  for (int j=2;j<=n;j++) se(j-1) = se(j);

  double f = 0.0;
  double tst1 = 0.0;
  se(n) = 0.0;
  for (int l=1;l<=n;l++)
  {
    int j = 0;
    double h = fabs(sd(l)) + fabs(se(l));
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
	double g = sd(l);
	double p = (sd(l1) - g) / (2.0 * se(l));
	double r = pythag(p,1.0);
	sd(l) = se(l) / (p + sign(r,p));
	sd(l1) = se(l) * (p + sign(r,p));
	double dl1 = sd(l1);
	h = g - sd(l);
	//		if (l2 <= n)
	for(int k=l2;k<=n;k++) sd(k) = sd(k) - h;
	f = f + h;
	//
	//     .......... ql transformation ..........
	//
	p = sd(m);
	double c = 1.0;
	double c2 = c, c3=c;
	double el1 = se(l1);
	double s = 0.0, s2=0.0;
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
    double p = sd(i);

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
