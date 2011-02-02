// File: tglrat.cc  Find eigenvalues of real symmetric tridiagonal matrix.

// Modifications for oml containers Copyright (1994-2003), Jan N. Reimers

#include "oml/vector.h"
#include <cassert>
#include <cmath>
#include <iostream>

double pythag(double,double);
double epsilon(double);

inline double sign(double a, double b) {return b<0 ? -fabs(a) : fabs(a);}

//--------------------------------------------------------------------------
//
//    this subroutine finds the eigenvalues of a symmetric
//    tridiagonal matrix by the rational ql method.
//
void tqlrat(Vector<double>& d, Vector<double>& e2 ,int ierr)
{
	int i,j,l,m,ii,l1,mml;
	double b=0,c=0,f,g,h,p,r,s,t;
	int n=d.size();
	
	std::cout << "in tqlrat" << std::endl;
	
	assert(d.size()==e2.size());

	ierr = 0;
	if (n == 1) return;

  for (i=2;i<=n;i++) e2(i-1) = e2(i);

	f = 0.0;
	t = 0.0;
	e2(n) = 0.0;

	for (l=1;l<=n;l++) //do 290 l = 1, n
	{
		j = 0;
		h = fabs(d(l)) + sqrt(e2(l));
		if (t <= h)
		{
			t = h;
			b = epsilon(t);
			c = b * b;
		}
//
//     .......... look for small squared sub-diagonal element ..........
//
    for (m=l;m<=n;m++) if (e2(m) <= c) break; //e2(n) is always zero, so there is no exit

    if (m != l)
		{
L130:
			if (j == 30) 
			{
				ierr=l;
				return;
			};
			j = j + 1;
			l1 = l + 1;
			s = sqrt(e2(l));
			g = d(l);
			p = (d(l1) - g) / (2.0 * s);
			r = pythag(p,1.0);
			d(l) = s / (p + sign(r,p));
			h = g - d(l);


			for (i=l1;i<=n;i++) d(i) = d(i) - h;
			f = f + h;
			g = d(m);
			if (g == 0.0) g = b;
			h = g;
			s = 0.0;
			mml = m - l;
			for (ii=1;ii<=mml;ii++)
			{
				i = m - ii;
				p = g * h;
				r = p + e2(i);
				e2(i+1) = s * r;
				s = e2(i) / r;
				d(i+1) = h + s * (h + d(i));
				g = d(i) - e2(i) / g;
				if (g == 0.0) g = b;
				h = g * p / r;
			}
			e2(l) = s * g;
			d(l) = h;

			if (h != 0.0 && fabs(e2(l)) > fabs(c/h)) //guard against underflow in convergence test
			{
				e2(l) = h * e2(l);
				if (e2(l) != 0.0) goto L130;
			}
		}
    p = d(l) + f;
//
//  order eigenvalues.
//
//		if (l == 1) {i=1;goto L270;}
		for (ii=2;ii<=l;ii++)
		{
			i = l + 2 - ii;
			if (p >= d(i-1)) goto L270;
			d(i) = d(i-1);
		}
    i = 1;
L270:
    d(i) = p;
  } //290 continue

	return;
}
