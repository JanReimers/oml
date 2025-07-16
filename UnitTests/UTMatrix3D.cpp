#include "gtest/gtest.h"
#include <iomanip>
#include <iostream>

using std::cout;
using std::endl;

import oml.Matrix3D;


class Matrix3DTests : public ::testing::Test
{
public:
    Matrix3DTests()
    {
        StreamableObject::SetToPretty();
    }
};

TEST_F(Matrix3DTests,Constructors)
{
    using mt=Matrix3D<double>;
    EXPECT_EQ(mt(),mt(1,0,0,0,1,0,0,0,1));
    mt M(0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9),M1;
    EXPECT_EQ(M1=M,mt(0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9));
}

TEST_F(Matrix3DTests,Operators)
{
    using mt=Matrix3D<double>;
    using vt=Vector3D<double>;
    mt m1(1.0,-1.0,0.5,-1,2,.5,1,.5,2),m2(0.5,1,1,2,1,.5,.5,2,2);
    vt v(1,2,2);
    Vector3D<int>    vi(1,2,2);

    EXPECT_EQ(m1+m2,mt(1.5,0,1.5,1,3,1,1.5,2.5,4 ));
    EXPECT_EQ(m1-m2,mt(0.5,-2,-0.5,-3,1,0,0.5,-1.5,0));
    EXPECT_EQ(m1+=m2,mt(1.5,0,1.5,1,3,1,1.5,2.5,4));
    EXPECT_EQ(m1-=m2,mt(1,-1,0.5,-1,2,0.5,1,0.5,2));
    EXPECT_EQ(m1*m2,mt(-1.25,1,1.5,3.75,2,1,2.5,5.5,5.25));
    EXPECT_EQ(m1*v,vt(0,4,6));
    EXPECT_EQ(v*m1,vt(1,4,5.5));
    EXPECT_EQ(v*m1*v,20);
    EXPECT_EQ(m1*vi,vt(0,4,6));
    EXPECT_EQ(vi*m1,vt(1,4,5.5));
    EXPECT_EQ(vi*m1*vi,20);
    EXPECT_EQ(m1*2.0,mt(2,-2,1,-2,4,1,2,1,4));
    EXPECT_EQ(2.0*m1,mt(2,-2,1,-2,4,1,2,1,4));
    EXPECT_EQ(m1/2.0,mt(0.5,-0.5,0.25,-0.5,1,0.25,0.5,0.25,1));
    EXPECT_EQ(~m1,mt(1,-1,0.5,-1,2,0.5,1,0.5,2));
    EXPECT_EQ(Determinant(m2),0.25);
    EXPECT_EQ(Invert(m2),mt(4,0,-2,-15,2,7,14,-2,-6));
    EXPECT_EQ(m2*Invert(m2),mt(1,0,0,0,1,0,0,0,1));
    EXPECT_TRUE(m1==m1);
    EXPECT_FALSE(m1==m2);
    EXPECT_FALSE(m1!=m1);
    EXPECT_TRUE(m1!=m2);
  }
