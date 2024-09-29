// File: UT-Vector3D-double.cc  Unit test the Vector3D class for double data types.

// Copyright (1994-2003), Jan N. Reimers

#include "gtest/gtest.h"
#include "oml/vector3d.h"
#include "oml/io3d.h"
#include "oml/imp/stream.h"
#include <cmath>
#include <iostream>
#include <sstream>
#include <fstream>
#include <iomanip>

using std::cout;
using std::endl;

class Vector3DTests : public ::testing::Test
{
public:
    Vector3DTests()
    {
        StreamableObject::SetToPretty();
    }
};

TEST_F(Vector3DTests,Constructors)
{
    EXPECT_EQ(Vector3D<double>().x,0.0);
    EXPECT_EQ(Vector3D<double>(0.1,0.2,0.3).x,0.1);
    EXPECT_EQ(Vector3D<double>(0.1,0.2,0.3).y,0.2);
    EXPECT_EQ(Vector3D<double>(0.1,0.2,0.3).z,0.3);
    Vector3D<double> V(0.3,0.2,0.1),V2;
    EXPECT_EQ(Vector3D<double>(V).x,0.3);
    EXPECT_EQ(Vector3D<double>(V).y,0.2);
    EXPECT_EQ(Vector3D<double>(V).z,0.1);
    V2=V;
    EXPECT_EQ(V2.x,0.3);
    EXPECT_EQ(V2.y,0.2);
    EXPECT_EQ(V2.z,0.1);
}

TEST_F(Vector3DTests,Operators)
{
//      cout.precision(16);
    typedef Vector3D<double> VT;
    Vector3D<double> v1(1.0,-1.0,0.5), v2(0.5,0.25,0.75);
    EXPECT_EQ(v1+v2,VT(1.5,-0.75,1.25));
    EXPECT_EQ(v1-v2,VT(0.5,-1.25,-0.25));
    EXPECT_EQ(v1+=v2,VT(1.5,-0.75,1.25));
    EXPECT_EQ(v1-=v2,VT(1,-1,0.5));
    EXPECT_EQ(v1*v2,0.625);
    EXPECT_EQ(v1*2.0,VT(2,-2,1));
    EXPECT_EQ(2.0*v1,VT(2,-2,1));
    EXPECT_EQ(v1/2.0,VT(0.5,-0.5,0.25));
    EXPECT_EQ(v1*=2.0,VT(2,-2,1));
    EXPECT_EQ(v1/=2.0,VT(1,-1,0.5));
    EXPECT_EQ(Cross(v1,v2),VT(-0.875,-0.5,0.75));
    EXPECT_EQ(Cross(v1,v2)*v1,0.);
    EXPECT_EQ(Cross(v1,v2)*v2,0.);
    EXPECT_EQ(norm(v1),1.5);
    EXPECT_EQ(normalize(v1),VT(0.6666666666666666,-0.6666666666666666,0.3333333333333333));
    EXPECT_EQ(+v1,VT(1,-1,0.5));
    EXPECT_EQ(-v1,VT(-1,1,-0.5));
    EXPECT_EQ(v1==v1,true);
    EXPECT_EQ(v1==v2,false);
    EXPECT_EQ(v1!=v1,false);
    EXPECT_EQ(v1!=v2,true);
    EXPECT_EQ(v1>v2,true);
    EXPECT_EQ(v1<v2,false);
    EXPECT_EQ(v1>=v2,true);
    EXPECT_EQ(v1>=v1,true);
    EXPECT_EQ(v1<=v2,false);
    EXPECT_EQ(v2<=v1,true);
    EXPECT_EQ(v1<=v1,true);
    EXPECT_EQ(angle(v1,v2),1.1091358114108441);
    EXPECT_EQ(angle_degrees(v1,v2),63.54880090065938);
    EXPECT_EQ(v1!=v2,true);
}

TEST_F(Vector3DTests,AsciiIO)
{
    StreamableObject::SetToAscii();
    Vector3D<double> A(.1,.2,.3),B;
    {
        std::ofstream file("temp.dat");
        file << A;
    }
    {
        std::ifstream file("temp.dat");
        file >> B;
    }
    StreamableObject::SetToPretty();
    EXPECT_EQ(B,A);
}
TEST_F(Vector3DTests,BinaryIO)
{
    StreamableObject::SetToBinary();
    Vector3D<double> A(.1,.2,.3),B;
    {
        std::ofstream file("temp.dat");
        file << A;
    }
    {
        std::ifstream file("temp.dat");
        file >> B;
    }
    StreamableObject::SetToPretty();
    EXPECT_EQ(B,A);
}


