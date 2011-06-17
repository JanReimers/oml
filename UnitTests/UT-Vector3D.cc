// File: UT-Vector3D-double.cc  Unit test the Vector3D class for double data types.

// Copyright (1994-2003), Jan N. Reimers

#include "oml/UnitTests/UnitTest.H"
#include "oml/vector3d.h"
#include "oml/io3d.h"
#include "oml/stream.h"
#include <cmath>
#include <iostream>
#include <strstream>
#include <fstream>
#include <iomanip>

using namespace std;

int TestVector3D()
{
  const char* Class="TVector3D<double>";
  bool pass=true;
//
//  Start message.
//

  StartClass(Class);
  StreamableObject::SetOutputMode(StreamableObject::pretty);

  {
    EXPECT(Vector3D<double>(),"(0,0,0)");
    EXPECT(Vector3D<double>(0.1,0.2,0.3),"(0.1,0.2,0.3)");
    Vector3D<double> V(0.3,0.2,0.1),V2;
    EXPECT(Vector3D<double>(V),"(0.3,0.2,0.1)");
    EXPECT3(V2,V2=V,"(0.3,0.2,0.1)");
  }
  {
    Vector3D<double> v1(1.0,-1.0,0.5), v2(0.5,0.25,0.75);
    EXPECT(v1,"(1,-1,0.5)");
    EXPECT(v2,"(0.5,0.25,0.75)");
    EXPECT(v1+v2,"(1.5,-0.75,1.25)");
    EXPECT(v1-v2,"(0.5,-1.25,-0.25)");
    EXPECT(v1+=v2,"(1.5,-0.75,1.25)");
    EXPECT(v1-=v2,"(1,-1,0.5)");
    EXPECT(v1*v2,"0.625");
    EXPECT(v1*2.0,"(2,-2,1)");
    EXPECT(2.0*v1,"(2,-2,1)");
    EXPECT(v1/2.0,"(0.5,-0.5,0.25)");
    EXPECT(v1*=2.0,"(2,-2,1)");
    EXPECT(v1/=2.0,"(1,-1,0.5)");
    EXPECT(Cross(v1,v2),"(-0.875,-0.5,0.75)");
    EXPECT(Cross(v1,v2)*v1,0.);
    EXPECT(Cross(v1,v2)*v2,0.);
    EXPECT(!v1,1.5);
    EXPECT(~v1,"(0.666667,-0.666667,0.333333)");
    EXPECT(+v1,"(1,-1,0.5)");
    EXPECT(-v1,"(-1,1,-0.5)");
    EXPECT(v1==v1,true);
    EXPECT(v1==v2,false);
    EXPECT(v1!=v1,false);
    EXPECT(v1!=v2,true);
    EXPECT(v1>v2,true);
    EXPECT(v1<v2,false);
    EXPECT(v1>=v2,true);
    EXPECT(v1>=v1,true);
    EXPECT(v1<=v2,false);
    EXPECT(v2<=v1,true);
    EXPECT(v1<=v1,true);
    EXPECT(v1|v2,"1.10914");
    EXPECT(v1||v2,"63.5488");
    EXPECT(v1!=v2,true);
  }

  StreamableObject::SetToAscii();
  {
    Vector3D<double> A(.1,.2,.3),B;
    EXPECT(A,"0.1 0.2 0.3 ");
    {
      ofstream file("temp.dat");
      file << A;
    }
    {
      ifstream file("temp.dat");
      file >> B;
    }
    StreamableObject::SetToPretty();
    EXPECT1(B,A,"file >> B");
  }
  StreamableObject::SetToBinary();
  {
    Vector3D<double> A(.1,.2,.3),B;
    {
      ofstream file("temp.dat");
      file << A;
    }
    {
      ifstream file("temp.dat");
      file >> B;
    }
    StreamableObject::SetToPretty();
    EXPECT1(B,A,"file >> B");
  }
  return pass ? 0 : -1;
}
