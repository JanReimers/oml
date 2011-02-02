// File: UT-SolveTriDiagonal.cc  Unit test the tri-diagonal system solver.

// Copyright (1994-2003), Jan N. Reimers

#include "oml/UnitTests/UnitTest.H"
#include "oml/numeric.h"
#include "oml/array.h"
#include "oml/vector.h"
#include "oml/vector_io.h"
#include "oml/matrix.h"
#include "oml/random.h"
#include <iostream>

int main()
{
  const char* Class="Tri diagonal system solver";
  StartClass(Class);
  StreamableObject::SetToPretty();
  int n=10;
  double a=0.2;
  bool pass=true;

  // Test std::vector version
  {
    VecLimits lim(0,n);
    // Load up diagonal and off diagonal elements.
    Vector<double> r(lim),off(lim),diag(lim);
    FillRandom(r,1.0);
    FillRandom(off,-a);
    FillRandomPositive(diag,1+2*a);
    // Load up a symmetrix tri diagonal matrix.
    Matrix<double> A(lim,lim);
    Fill(A,0.0);
    A.GetDiagonal()=diag;    
    for (int i=lim.Low;i<=lim.High-1;i++) 
      {
	A(i,i+1)=off(i);
	A(i+1,i)=off(i);
      }
    
    std::cout.setf(std::ios::fixed,std::ios::floatfield);
    std::cout << "A=" << A << std::endl;
    
    Vector<double> u=SolveTriDiag(off,diag,r);
    std::cout << "A*u-r=" << A*u-r << std::endl;
    EXPECT(!(A*u-r) < 1e-6,"1");
  }

  // Test array version.  This on is ugly because we can't multiply an Array * Matrix.
  {
    // Load up diagonal and off diagonal elements.
    Array<double> r(n),off(n),diag(n);
    FillRandom(r,1.0);
    FillRandom(off,a);
    FillRandomPositive(diag,1+2*a);
    std::cout.setf(std::ios::fixed,std::ios::floatfield);
    
    Array<double> u=SolveTriDiag(off,diag,r);

    double err=0.0;
    std::cout << "T*u-r=";
    for (int i=0;i<n;i++)
      {
	double t=u[i]*diag[i]-r[i];
	if (i>0) t+=u[i-1]*off[i-1];
	if (i<n-1) t+=u[i+1]*off[i];
	std::cout << t << " ";
	err+=t*t;
      }
    std::cout << std::endl;
    EXPECT(sqrt(err) < 1e-6,"1");
  }

  return pass ? 0 : -1;
}
