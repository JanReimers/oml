// File: UT-Matrix3D-double.cc  Unit test the Matrix3D class for double data types.

// Copyright (1994-2005), Jan N. Reimers

#include "oml/UnitTests/UnitTest.H"
#include "oml/matrix3d.h"
#include "oml/io3d.h"
#include "oml/stream.h"
#include <iostream>
#include <iomanip>
#include <fstream>

int TestMatrix3DDouble()
{
  const char* Class="Matrix3D<double>";
  bool pass=true;
//
//  Start message.
//

  StartClass(Class);
  StreamableObject::SetOutputMode(StreamableObject::pretty);
#if __GNUG__>=3
  std::cout.setf(std::_Ios_Fmtflags(0),std::ios::floatfield);
#else
  std::cout.setf(0,std::ios::floatfield);
#endif
  std::cout.precision(6);

//###########################################################################
//
//  Constructors
//
  {
    EXPECT(Matrix3D<double>(),"[ 1 0 0 ]\n[ 0 1 0 ]\n[ 0 0 1 ]\n");
    EXPECT(Matrix3D<double>(0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9),"[ 0.1 0.2 0.3 ]\n[ 0.4 0.5 0.6 ]\n[ 0.7 0.8 0.9 ]\n");
    Matrix3D<double> M(0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9),M1;
    EXPECT(Matrix3D<double>(M),"[ 0.1 0.2 0.3 ]\n[ 0.4 0.5 0.6 ]\n[ 0.7 0.8 0.9 ]\n");
    EXPECT(M1=M,"[ 0.1 0.2 0.3 ]\n[ 0.4 0.5 0.6 ]\n[ 0.7 0.8 0.9 ]\n");
  }
  {
    Matrix3D<double> m1(1.0,-1.0,0.5,  -1,2,.5,  1,.5,2), m2(0.5,1,1,  2,1,.5,  .5,2,2);
    Vector3D<double> v(1,2,2);
    Vector3D<int>    vi(1,2,2);

    EXPECT(m1,"[ 1 -1 0.5 ]\n[ -1 2 0.5 ]\n[ 1 0.5 2 ]\n");
    EXPECT(m1+m2,"[ 1.5 0 1.5 ]\n[ 1 3 1 ]\n[ 1.5 2.5 4 ]\n");
    EXPECT(m1-m2,"[ 0.5 -2 -0.5 ]\n[ -3 1 0 ]\n[ 0.5 -1.5 0 ]\n");
    EXPECT(m1+=m2,"[ 1.5 0 1.5 ]\n[ 1 3 1 ]\n[ 1.5 2.5 4 ]\n");
    EXPECT(m1-=m2,"[ 1 -1 0.5 ]\n[ -1 2 0.5 ]\n[ 1 0.5 2 ]\n");
    EXPECT(m1*m2,"[ -1.25 1 1.5 ]\n[ 3.75 2 1 ]\n[ 2.5 5.5 5.25 ]\n");
    EXPECT(m1*v,"(0,4,6)");
    EXPECT(v*m1,"(1,4,5.5)");
    EXPECT(v*m1*v,20);
    EXPECT(m1*vi,"(0,4,6)");
    EXPECT(vi*m1,"(1,4,5.5)");
    EXPECT(vi*m1*vi,20);
    EXPECT(m1*2.0,"[ 2 -2 1 ]\n[ -2 4 1 ]\n[ 2 1 4 ]\n");
    EXPECT(2.0*m1,"[ 2 -2 1 ]\n[ -2 4 1 ]\n[ 2 1 4 ]\n");
    EXPECT(m1/2.0,"[ 0.5 -0.5 0.25 ]\n[ -0.5 1 0.25 ]\n[ 0.5 0.25 1 ]\n");
    EXPECT(~m1,"[ 1 -1 0.5 ]\n[ -1 2 0.5 ]\n[ 1 0.5 2 ]\n");
    EXPECT(Determinant(m2),0.25);
    EXPECT(Invert(m2),"[ 4 0 -2 ]\n[ -15 2 7 ]\n[ 14 -2 -6 ]\n");
    EXPECT(m2*Invert(m2),"[ 1 0 0 ]\n[ 0 1 0 ]\n[ 0 0 1 ]\n");
    EXPECT(m1==m1,true);
    EXPECT(m1==m2,false);
    EXPECT(m1!=m1,false);
    EXPECT(m1!=m2,true);
  }

  StreamableObject::SetToAscii();
  {
    Matrix3D<double> A(.1,.2,.3,.4,.5,.6,.7,.8,.9),B;
    {
      std::ofstream file("temp.dat");
      file << A;
    }
    {
      std::ifstream file("temp.dat");
      file >> B;
    }
      StreamableObject::SetToPretty();
      EXPECT1(B,A,"file >> B");
  }
  StreamableObject::SetToBinary();
  {
    Matrix3D<double> A(.1,.2,.3,.4,.5,.6,.7,.8,.9),B;
    {
      std::ofstream file("temp.dat");
      file << A;
    }
    {
      std::ifstream file("temp.dat");
      file >> B;
    }
      StreamableObject::SetToPretty();
      EXPECT1(B,A,"file >> B");
  }
  return pass ? 0 : -1;
}
