// File: UT-List-complex.cc  Unit test the List class for complex data types.

// Copyright (1994-2005), Jan N. Reimers

#include "oml/UnitTests/UnitTest.H"
#include "oml/list.h"
#include "oml/list_io.h"
#include "oml/array_io.h"
#include "oml/omlcomplex.h"
#include <fstream>


int main()
{
  const char* Class="List<std::complex<double> >";
  bool pass=true;
//
//  Start message.
//

  StartClass(Class);
  StreamableObject::SetOutputMode(StreamableObject::pretty);
  {
    List<std::complex<double> > A,B;
    std::complex<double> a(1,M_PI);
    std::complex<double> b(0.25,M_PI/4);
    for (int i=1;i<10;i++) A.push_back(a);
    for (int i=1;i<10;i++) B.push_back(b);  
    
    EXPECT(Max(abs(std::complex<double> (1.0)/A   -(std::complex<double>(1.0)/ a)))<eps,true);
    EXPECT(Max(abs(sin (A) - sin(a)))<eps,true);
    EXPECT(Max(abs(cos (A) - cos(a)))<eps,true);
    EXPECT(Max(abs(sinh(A) -sinh(a)))<eps,true);
    EXPECT(Max(abs(cosh(A) -cosh(a)))<eps,true);
    EXPECT(Max(abs(log (A) -log (a)))<eps,true);
    EXPECT(Max(abs(exp (A) -exp (a)))<eps,true);
    EXPECT(Max(abs(sqrt(A) -sqrt(a)))<eps,true);
    EXPECT(Max(abs(pow2(A) -     a*a))<10 *eps,true);
    EXPECT(Max(abs(pow3(A) -   a*a*a))<10 *eps,true);
    EXPECT(Max(abs(pow4(A) - a*a*a*a))<100*eps,true);
    EXPECT(Max(abs(pow (A,std::complex<double>(0.6))-pow   (a,0.6)))<4*eps,true);
    EXPECT(Max(abs(pow (A,B)  -pow  (a,b)))<eps,true);
    EXPECT(Max(abs(conj(A) -conj(a)))<eps,true);
    EXPECT(Max(abs(norm(A) -norm(a)))<2*eps,true);
    EXPECT(Max(abs(arg (A) -arg(a)))<eps,true);
    EXPECT(Max(abs(real(A) -real(a)))<eps,true);
    EXPECT(Max(abs(imag(A) -imag(a)))<eps,true);
  }
	
//#########################################################################
//
//  Test data file IO.
//
  {
    StreamableObject::SetToAscii();
    List<std::complex<double> > A,B;
    for (int i=1;i<10;i++) A.push_back(i);
    {
      std::ofstream file("temp.dat");
      file << A;
    }

    {
      std::ifstream file("temp.dat");
      //    EXPECT3(B,file >> B,A);
    }

    StreamableObject::SetToBinary();
    {
      std::ofstream file("temp.dat");
      file << A;
    }

    {
      std::ifstream file("temp.dat");
      //    EXPECT3(B,file >> B,A);
    }
  }
  return pass ? 0 : -1;
}


