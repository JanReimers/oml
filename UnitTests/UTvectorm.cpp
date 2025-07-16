#include "gtest/gtest.h"
#include <iomanip>
#include <iostream>
#include <fstream>
#include <complex>

import oml.Vector;

template <class T> class VectorTests : public ::testing::Test
{
public:
    VectorTests()
    {
        StreamableObject::SetToPretty();
    }
};

class VectorRealTests : public ::testing::Test
{
public:
    VectorRealTests()
    {
        StreamableObject::SetToPretty();
    }
};

class VectorComplexTests : public ::testing::Test
{
public:
    VectorComplexTests()
    {
        StreamableObject::SetToPretty();
    }
};

TYPED_TEST_SUITE_P(VectorTests);

TYPED_TEST_P(VectorTests,Constructors)
{
    typedef Vector<TypeParam> VectorT;
    VectorT V0;
    EXPECT_EQ(V0.size(),0);
    VectorT V1=VectorT();
    EXPECT_EQ(V1.size(),0);
}


TYPED_TEST_P(VectorTests,OverloadedOperators1)
{
    typedef Vector<TypeParam> VectorT;
    int N=10;
    VectorT V1(N),V2(N);
    Fill<TypeParam>(V1,1.0);
    Fill<TypeParam>(V2,1.0);
    EXPECT_EQ(V1,V2);
    EXPECT_EQ(V1==V2,true);
    EXPECT_EQ(V1!=V2,false);
    EXPECT_EQ((V1+V2)(1),2.0); // Use Fill constructors EXPECT_EQ((A+B),MatrixT(A.GetLimits(),2.0)
    //EXPECT_EQ((A+B),MatrixT(A.GetLimits(),2.0); TODO
    EXPECT_EQ((V1-V2  )(1),0.0);
    EXPECT_EQ((V1*3.0 )(1),3.0);
    EXPECT_EQ((V1*3   )(1),3.0);
    EXPECT_EQ((3.0*V1 )(1),3.0);
    EXPECT_EQ((3*V1   )(1),3.0);
    EXPECT_EQ((V1/2.0 )(1),0.5);
    EXPECT_EQ((V1/2   )(1),0.5);
    EXPECT_EQ((V1+=V2 )(1),2.0);
    EXPECT_EQ((V1-=V2 )(1),1.0);
    EXPECT_EQ((V1*=3.0)(1),3.0);
    EXPECT_EQ((V1*=3  )(1),9.0);
    EXPECT_EQ((V1/=3.0)(1),3.0);
    EXPECT_EQ((V1/=3  )(1),1.0);

    EXPECT_EQ((V1+V2-V1*2+V2*V1)(1),10.0);

}

TYPED_TEST_P(VectorTests,AsciiAndBinaryIO)
{
    typedef Vector<TypeParam> VectorT;
    StreamableObject::SetToAscii();
    {
        VectorT A(VecLimits(-10,20)),B;
        FillRandom(A,TypeParam(100));
        {
            std::ofstream file("temp.dat");
            file << std::setprecision(12) << A;
        }
        StreamableObject::SetToPretty();
        {
            std::ifstream file("temp.dat");
            file >> B;
        }
        EXPECT_NEAR(Max(fabs(B-A)),0.0,1e-10);
    }

    StreamableObject::SetToBinary();
    {
        VectorT A(VecLimits(-10,20)),B;
        FillRandom(A,TypeParam(100));
        {
            std::ofstream file("tempb.dat",std::ios::out|std::ios::binary);
            file << A;
        }
        StreamableObject::SetToPretty();
        {
            std::ifstream file("tempb.dat",std::ios::in|std::ios::binary);
            file >> B;
        }
        EXPECT_EQ(B,A);
    }
}

TEST_F(VectorRealTests,MinMax)
{
    typedef Vector<double> VectorT;
    int N=8;
    VectorT A(N),B(A);
    FillLinear(A,1.0/4,2.0); //all elements should be exactly represented in floating point
    FillLinear(B,1.0/4,2.0); //all elements should be exactly represented in floating point

    EXPECT_EQ(Max(A),2.0);
    EXPECT_EQ(Min(A),0.25);
    EXPECT_EQ(Max(A+B),4.0);
    EXPECT_EQ(Min(A+B),0.5);

    {
        int N=1000;
        VectorT Ran(N);
        FillRandom(Ran,100.0);
        double min=Ran(1), max=Ran(1);
        for (int i=1; i<=N; i++)
            {
                min=Min(min,Ran(i));
                max=Max(max,Ran(i));
            }
        EXPECT_EQ(Min(Ran),min);
        EXPECT_EQ(Max(Ran),max);
    }

}

TEST_F(VectorComplexTests,fabsConj)
{
    typedef std::complex<double> dcmplx;
    typedef Vector <dcmplx> VectorCT;

    int N=10;
    VectorCT A(N),B(A);
    Fill(A,dcmplx(2,-2));
    B=conj(A);
    EXPECT_EQ(real(Sum(A+B)),40.);
    EXPECT_EQ(Sum(real(A+B)),40.);
    EXPECT_EQ(imag(Sum(A+B)),  0.);
    EXPECT_EQ(Sum(imag(A+B)),  0.);
    EXPECT_EQ(real(A*B),80.);
    EXPECT_EQ(imag(A*B), 0.);

}

TEST_F(VectorComplexTests,MixedTypes)
{
    typedef std::complex<double> dcmplx;
    typedef Vector <dcmplx> VectorCT;
    typedef Vector <double> VectorRT;

    int N=10;
    VectorCT Vc(N,dcmplx(0.25,0.5));
    VectorRT Vr(VecLimits(N),2.0);
    dcmplx   Sc(0.5,-0.25);
    double   Sr(4.0);
    // Try every combination we can think of
    EXPECT_EQ((Vc*Vc)     ,dcmplx(-1.875,2.5));
    EXPECT_EQ((Vc*Vr)     ,dcmplx(5,10));
    EXPECT_EQ((Vr*Vc)     ,dcmplx(5,10));
    EXPECT_EQ((Vc*Sr)(1)  ,dcmplx(1,2));
    EXPECT_EQ((Sr*Vc)(1)  ,dcmplx(1,2));
    EXPECT_EQ((Sc*Vr)(1)  ,dcmplx(1,-0.5));
    EXPECT_EQ((Vr*Sc)(1)  ,dcmplx(1,-0.5));

    VectorRT V1=-0.5*imag(Vc-Vc);
    EXPECT_EQ(V1(1),0.0);
}


template <class T, typename Tf> void TestUnopLight(T dummy,Tf f)
{
//    cout << __PRETTY_FUNCTION__ << endl;
    typedef Vector <T> VectorT;
    index_t N=10;
    T s(0.5);
    VectorT V(VecLimits(N),s);
    Fill(V,s);
    EXPECT_EQ(f(V),f(s));
}

template <class T, typename Tf> void TestUnop(T dummy,Tf f)
{
//    cout << __PRETTY_FUNCTION__ << endl;
    typedef Vector <T> VectorT;
    index_t N=10;
    T s(0.5);
    VectorT V(VecLimits(N),s);
    EXPECT_EQ(f(V)(1  ),f(s));
//    EXPECT_NEAR(fabs(f(V+V)(1)-f(2.*s)),0.0,eps);
}

TYPED_TEST_P(VectorTests,UnaryOps)
{
    TestUnop(TypeParam(0),[](const auto &x) { return sin(x); });
    TestUnop(TypeParam(0),[](const auto &x) { return cos(x); });
    TestUnop(TypeParam(0),[](const auto &x) { return tan(x); });
    //TestUnop(TypeParam(0),[](const auto &x) { return asin(x); }); Restricted range make these tricky
    //TestUnop(TypeParam(0),[](const auto &x) { return acos(x); }); Uncomment to make sure they at least compile
    //TestUnop(TypeParam(0),[](const auto &x) { return atan(x); });
    TestUnop(TypeParam(0),[](const auto &x) { return sinh(x); });
    TestUnop(TypeParam(0),[](const auto &x) { return cosh(x); });
    TestUnop(TypeParam(0),[](const auto &x) { return tanh(x); });
    TestUnop(TypeParam(0),[](const auto &x) { return exp(x); });
    TestUnop(TypeParam(0),[](const auto &x) { return log(x); });
    TestUnop(TypeParam(0),[](const auto &x) { return log10(x); });
    //TestUnop(TypeParam(0),[](const auto &x) { return pow10(x); });
    TestUnop(TypeParam(0),[](const auto &x) { return sqrt(x); });
    TestUnop(TypeParam(0),[](const auto &x) { return fabs(x); });
    TestUnop(TypeParam(0),[](const auto &x) { return -x; });
    TestUnop(TypeParam(0),[](const auto &x) { return +x; });

    TestUnopLight(TypeParam(0),[](const auto &x) { return isnan(x); });
    TestUnopLight(TypeParam(0),[](const auto &x) { return isinf(x); });

}

TEST_F(VectorComplexTests,UnaryOps)
{
   typedef std::complex<double> dcmplx;
   TestUnop(dcmplx(0),[](const auto &x) { return conj(x); });
   TestUnop(dcmplx(0),[](const auto &x) { return real(x); });
   TestUnop(dcmplx(0),[](const auto &x) { return imag(x); });
   TestUnop(dcmplx(0),[](const auto &x) { return norm(x); });
   TestUnop(dcmplx(0),[](const auto &x) { return arg(x); });

}

TYPED_TEST_P(VectorTests,BinaryOps)
{
    typedef Vector <TypeParam> VectorT;
    index_t N=10;
    TypeParam sA(0.5),sB(2.0);
    VectorT A(VecLimits(N),sA),B(VecLimits(N),sB);
    EXPECT_EQ((A+B)(1),sA+sB);
    EXPECT_EQ((A-B)(1),sA-sB);
    EXPECT_EQ(DirectMultiply(A,B)(1),sA*sB);
    EXPECT_EQ(DirectDivide  (A,B)(1),sA/sB);
    EXPECT_EQ((A+sB)(1),sA+sB);
    EXPECT_EQ((A-sB)(1),sA-sB);
    EXPECT_EQ((sA+B)(1),sA+sB);
    EXPECT_EQ((sA-B)(1),sA-sB);
    EXPECT_EQ((A*sB)(1),sA*sB);
    EXPECT_EQ((A/sB)(1),sA/sB);
    EXPECT_EQ((sA*B)(1),sA*sB);
    EXPECT_TRUE(A==sA);
    EXPECT_TRUE(sA==A);
    EXPECT_TRUE(A==A);
    EXPECT_TRUE(A!=B);
    EXPECT_TRUE(A!=sB);
    EXPECT_TRUE(sB!=A);

}


TEST_F(VectorComplexTests,RangeBasedLoops)
{
    using dcmplx=std::complex<double>;
    typedef Vector<dcmplx> VectorCT;
    typedef Vector<double> VectorRT;

    int N=10;
    VectorCT Vc(VecLimits(N),dcmplx(0.25,0.5));
    VectorRT Vr(VecLimits(N),2.0);
    dcmplx   Sc(0.5,-0.25);
    double   Sr(4.0);

    for (double d:Vr) Sr+=d;
    for (dcmplx c:Vc) Sc+=c;

}


REGISTER_TYPED_TEST_SUITE_P(
            VectorTests,
            Constructors,
            OverloadedOperators1,
            AsciiAndBinaryIO,
            UnaryOps,
            BinaryOps
            );
//using MyTypes = ::testing::Types<double>;
using MyTypes = ::testing::Types<double,std::complex<double>>;
INSTANTIATE_TYPED_TEST_SUITE_P(My, VectorTests, MyTypes);
