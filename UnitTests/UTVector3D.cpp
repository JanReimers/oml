#include "gtest/gtest.h"
#include <iomanip>
import oml.Vector3D;
 

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
    EXPECT_EQ(Vector3D<double>(),Vector3D<double>(0,0,0));
    Vector3D<double> V(0.3,0.2,0.1),V2;
    EXPECT_EQ(V2=V,Vector3D<double>(0.3,0.2,0.1));
}


TEST_F(Vector3DTests,Operators)
{
    using vt=Vector3D<double>;
    vt v1(1.0,-1.0,0.5), v2(0.5,0.25,0.75);
    EXPECT_EQ(v1+v2,vt(1.5,-0.75,1.25));
    EXPECT_EQ(v1-v2,vt(0.5,-1.25,-0.25));
    EXPECT_EQ(v1+=v2,vt(1.5,-0.75,1.25));
    EXPECT_EQ(v1-=v2,vt(1,-1,0.5));
    EXPECT_EQ(v1*v2,0.625);
    EXPECT_EQ(v1*2.0,vt(2,-2,1));
    EXPECT_EQ(2.0*v1,vt(2,-2,1));
    EXPECT_EQ(v1/2.0,vt(0.5,-0.5,0.25));
    EXPECT_EQ(v1*=2.0,vt(2,-2,1));
    EXPECT_EQ(v1/=2.0,vt(1,-1,0.5));
    EXPECT_EQ(Cross(v1,v2),vt(-0.875,-0.5,0.75));
    EXPECT_EQ(Cross(v1,v2)*v1,0.);
    EXPECT_EQ(Cross(v1,v2)*v2,0.);
    EXPECT_EQ(norm(v1),1.5);
    EXPECT_EQ(normalize(v1),vt(2./3.,-2./3.,1./3.));
    EXPECT_EQ(+v1,vt(1,-1,0.5));
    EXPECT_EQ(-v1,vt(-1,1,-0.5));
    EXPECT_TRUE (v1==v1);
    EXPECT_FALSE(v1==v2);
    EXPECT_FALSE(v1!=v1);
    EXPECT_TRUE(v1!=v2);
    EXPECT_TRUE(v1>v2);
    EXPECT_FALSE(v1<v2);
    EXPECT_TRUE(v1>=v2);
    EXPECT_TRUE(v1>=v1);
    EXPECT_FALSE(v1<=v2);
    EXPECT_TRUE(v2<=v1);
    EXPECT_TRUE(v1<=v1);
    EXPECT_EQ(angle(v1,v2),1.1091358114108441);
    EXPECT_EQ(angle_degrees(v1,v2),63.54880090065938);
    EXPECT_TRUE(v1!=v2);
  }

