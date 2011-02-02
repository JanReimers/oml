// File: UT-Array-complex.cc  Unit test the Array class for complex data types.

// Copyright (1994-2005), Jan N. Reimers

#include "oml/UnitTests/UnitTest.H"
#include "oml/array.h"
#include "oml/array_io.h"
#include "oml/omlcomplex.h"
#include "oml/random.h"
#include <cmath>
#include <fstream>

int main()
{
  const char* Class="Array<std::complex<T> >";
  bool pass=true;
//
//  Start message.
//

  StartClass(Class);
  StreamableObject::SetToPretty();

  Array<std::complex<double> > A(10),B(10);
  std::complex<double> a(1,M_PI);
  std::complex<double> b(0.25,M_PI/4);
  Fill(A,a);
  Fill(B,b);

  EXPECT(Max(abs(std::complex<double> (1.0)/A   -(std::complex<double>(1.0)/ a)))<eps,true);
  EXPECT(Max(abs((A*=2.0)   -a*2.0))<eps,true);
  EXPECT(Max(abs((A/=2.0)   -a    ))<eps,true);
  EXPECT(Max(abs(sin(A)  - sin (a)))<eps,true);
  EXPECT(Max(abs(cos(A)  - cos (a)))<eps,true);
  EXPECT(Max(abs(sinh(A) -sinh (a)))<eps,true);
  EXPECT(Max(abs(cosh(A) -cosh (a)))<eps,true);
  EXPECT(Max(abs(log(A)  -log  (a)))<eps,true);
  EXPECT(Max(abs(exp(A)  -exp  (a)))<eps,true);
  EXPECT(Max(abs(sqrt(A) -sqrt (a)))<eps,true);
  EXPECT(Max(abs(pow2(A) -    a*a))<10*eps,true);
  EXPECT(Max(abs(pow3(A) - a*a*a))<10*eps,true);
  EXPECT(Max(abs(pow4(A)  -a*a*a*a))<100*eps,true);
  EXPECT(Max(abs(pow(A,std::complex<double>(0.6))-pow   (a,0.6)))<4*eps,true);
  EXPECT(Max(abs(pow(A,B)  -pow  (a,b)))<eps,true);
  EXPECT(Max(abs(conj(A) -conj(a)))<eps,true);
  EXPECT(Max(abs(norm(A) -norm(a)))<2*eps,true);
  EXPECT(Max(abs(arg(A) -arg(a)))<eps,true);
  EXPECT(Max(abs(real(A) -real(a)))<eps,true);
  EXPECT(Max(abs(imag(A) -imag(a)))<eps,true);
	
//-----------------------------------------------------
//
//  IO tests
//	
  StreamableObject::SetToAscii();
  {
    Array<std::complex<double> > A(1000),B;
    FillRandom(A,std::complex<double>(100.));
    {
      std::ofstream file("temp.dat");
      file << A;
    }
    StreamableObject::SetToPretty();
    {
      std::ifstream file("temp.dat");
      file >> B;
    }
    EXPECT1(B,A,"file >> B in ascii mode");
  }
  
  StreamableObject::SetToBinary();
  {
    Array<std::complex<double> > A(1000),B;
    FillRandom(A,std::complex<double>(100.));
    {
      std::ofstream file("tempb.dat",std::ios::out|std::ios::binary);
      file << A;
    }
    StreamableObject::SetToPretty();
    {
      std::ifstream file("tempb.dat",std::ios::in|std::ios::binary);
      file >> B;
    }
    EXPECT1(B,A,"file >> B in binary mode");
  }		
  return pass ? 0 : -1;
}

