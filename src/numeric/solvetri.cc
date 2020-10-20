// File TriDiag.C  Solve a tri diagonal cyctem arising from the implicit diff eqn.

#include "oml/vector.h"
#include <iostream>
#include <cstdlib>

// Copyright (1994-2003), Jan N. Reimers

//---------------------------------------------------------------------
//
//  Solve [Tri][u]=[r] for u.  Where Tri is symmetric and looks like this
//
//   b1 a1 0  0  0  0 0 . 0    0    0
//   a1 b2 a2 0  0  0 0 . 0    0    0
//   0  a2 b3 a3 0  0 0 . 0    0    0
//   0  0  a3 b4 a4 0 0 . 0    0    0
//   . . . . . . . . . . .
//   . . . . . . . . . . .
//   0  0  0  0  0  0 0   aN-2 bN-1 aN-1
//   0  0  0  0  0  0 0 . 0    aN-1 bN
//
// Note that aN is never used.
// See Numerical recipes section 2.4 (non-symmetric case).
//
Array<double> SolveTriDiag(const Array<double> & OffDiagonal,
			   const Array<double> & Diagonal,
			   const Array<double> & r)
{
  assert(OffDiagonal.size()==Diagonal.size());
  assert(OffDiagonal.size()==r.size());
  assert(Diagonal[0]!=0.0);

  int n=OffDiagonal.size();
  Array<double> gam(n),u(n); //Work space and return array.
  double beta=Diagonal[0];
  u[0]=r[0]/beta;
  for (int j=1;j<n;j++)
  {
    gam[j]=OffDiagonal[j-1]/beta;
    beta=Diagonal[j]-OffDiagonal[j-1]*gam[j];
    if (beta==0.0)
    {
      std::cerr << "TriDiag failed j=" << j << std::endl;
      exit(-1);
    }
    u[j]=(r[j]-OffDiagonal[j-1]*u[j-1])/beta;
  }
  for (int j=n-2;j>=0;j--)
    u[j]-=gam[j+1]*u[j+1];

  return u;
}

Vector<double> SolveTriDiag(const Vector<double> & OffDiagonal,
			    const Vector<double> & Diagonal,
			    const Vector<double> & r)
{
  assert(OffDiagonal.GetLimits()==Diagonal.GetLimits());
  assert(OffDiagonal.GetLimits()==r       .GetLimits());
  int low=r.GetLow();
  int high=r.GetHigh();

  assert(Diagonal(low)!=0.0);

  VecLimits l=OffDiagonal.GetLimits();
  Vector<double> gam(l),u(l); //Work space.
  double beta=Diagonal(low);
  u(low)=r(low)/beta;
  for (int j=low+1;j<=high;j++)
  {
    gam(j)=OffDiagonal(j-1)/beta;
    beta=Diagonal(j)-OffDiagonal(j-1)*gam(j);
    if (beta==0.0)
    {
      std::cerr << "TriDiag failed j=" << j << std::endl;
      exit(-1);
    }
    u(j)=(r(j)-OffDiagonal(j-1)*u(j-1))/beta;
  }
  for (int j=high-1;j>=low;j--)
    u(j)-=gam(j+1)*u(j+1);

  return u;
}
