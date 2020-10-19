// File: UT-DMatrix-double.cc  Unit test the DMatrix class for double data types.

// Copyright (1994-2020), Jan N. Reimers

#include "oml/dmatrix.h"
#include "oml/random.h"
#include "oml/stopw.h"
#include "gtest/gtest.h"
#include <iostream>
#include <fstream>
#include <complex>

using std::cout;
using std::endl;

//#include "oml/vector_io.h"
//#include "oml/matrix.h"
//#include "oml/minmax.h"


template <class T> class MatrixTests : public ::testing::Test
{
public:
    MatrixTests()
    {
        StreamableObject::SetToPretty();
    }
};

class MatrixRealTests : public ::testing::Test
{
public:
    MatrixRealTests()
    {
        StreamableObject::SetToPretty();
    }
};

class MatrixComplexTests : public ::testing::Test
{
public:
    MatrixComplexTests()
    {
        StreamableObject::SetToPretty();
    }
};

template <class T> void mmul(     DMatrix<T>& C,const DMatrix<T>& A, const DMatrix<T>& B);
template <class T> T    vmul(const Vector<T> V1,const DMatrix<T>& A, const Vector<T>& V2);

//
//  Hand coded
//
TYPED_TEST_SUITE_P(MatrixTests);

TYPED_TEST_P(MatrixTests,Constructors)
{
    typedef DMatrix<TypeParam> MatrixT;
    MatrixT A0;
    EXPECT_EQ(A0.GetNumRows(),0);
    EXPECT_EQ(A0.GetNumCols(),0);

    MatrixT A1(10,5);
    EXPECT_EQ(A1.GetNumRows(),10);
    EXPECT_EQ(A1.GetNumCols(),5);

    MatLimits lim(0,5,0,10);
    MatrixT A2(lim);
    EXPECT_EQ(A2.GetNumRows(),6);
    EXPECT_EQ(A2.GetNumCols(),11);
    EXPECT_EQ(A2.GetLimits(),lim);
    EXPECT_EQ(A2.GetRowLow(),0);
    EXPECT_EQ(A2.GetColLow(),0);
    EXPECT_EQ(A2.size(),66);

    MatrixT A3(A2),A4,A5(5,4);
    EXPECT_EQ(A3.GetLimits(),A2.GetLimits());
    EXPECT_NE(A4.GetLimits(),A2.GetLimits());
    EXPECT_NE(A5.GetLimits(),A2.GetLimits());
    A4=A2;
    A5=A2;
    EXPECT_EQ(A4.GetLimits(),A2.GetLimits());
    EXPECT_EQ(A5.GetLimits(),A2.GetLimits());
}

TYPED_TEST_P(MatrixTests,Fill_SetLimits_SubMatrix)
{
    typedef DMatrix<TypeParam> MatrixT;
    MatrixT A1(10,5);
    FillRandom(A1);
    MatrixT A2=A1.SubMatrix(MatLimits(5,4));
    A1.SetLimits(5,4,true);
    EXPECT_EQ(A1,A2);

    typedef Vector <TypeParam> VectorT;
    VectorT V1(10,TypeParam(0.5));
 //   VectorT V1(10,0.5); for T=complex this fails since the compiler promotes 0.5 to int rather than complex.
// using  enum  FillType; C++20
    VectorT V2(10,FillType::Zero);
    VectorT V3(10,FillType::Random);
    VectorT V4(10,FillType::Unit);
    VectorT V5(10);
    EXPECT_EQ(V1(10),0.5);
    EXPECT_EQ(V2(10),0.0);
    EXPECT_EQ(V4(10),1.0);

    V5.Fill(2.0);
    EXPECT_EQ(V5(10),2.0);
    V5.FillRandom();
}

TYPED_TEST_P(MatrixTests,Transpose_Slices)
{
    typedef DMatrix<TypeParam> MatrixT;
    typedef Vector <TypeParam> VectorT;
    MatrixT A1(10,5);
    FillRandom(A1);
    MatrixT A2=Transpose(A1);
    EXPECT_EQ(A2.GetLimits(),MatLimits(5,10));
    EXPECT_EQ(A1(10,5),A2(5,10));
    EXPECT_EQ(A1.GetRow(2),A2.GetColumn(2));
//    A1.Transpose(); TODO: transpose in place
//    EXPECT_EQ(A1.GetLimits(),A2.GetLimits());

    index_t N=10;
    MatrixT A3(N,N);
    FillRandom(A3);
    VectorT D=A3.GetDiagonal();
    EXPECT_EQ(D.size(),N);
    for (index_t i=1; i<=N; i++)
        EXPECT_EQ(A3(i,i),D(i));
}

TYPED_TEST_P(MatrixTests,SumDotMaxMinDirectMultiply)
{
    typedef DMatrix<TypeParam> MatrixT;
    index_t M=2,N=4,mn=M*N;
    MatrixT A(M,N),B(A);
    FillLinear<TypeParam>(A,1.0/4,2.0); //all elements should be exactly represented in floating point
    //cout << A << endl;
    EXPECT_EQ(Sum(A),mn*(mn+1)/2/4.);
    EXPECT_EQ(Dot(A,A),51.0/4); //Look up sum n^2 rule
    Fill<TypeParam>(B,0.5);
    EXPECT_EQ(DirectMultiply(A,B)(1,1),1.0/8);
    EXPECT_EQ(DirectDivide(A,B)(1,1),1.0/2);

    A.SetLimits(N,N);
    Unit(A);
    EXPECT_EQ(A,~A);
    EXPECT_EQ(A(2,2),1.0);
    EXPECT_EQ(A(2,1),0.0);
    EXPECT_EQ(Sum(A),TypeParam(N));
    EXPECT_EQ(IsSymmetric(A),true);
    A(1,2)=-1.0;
    EXPECT_EQ(IsSymmetric(A),false);
    MakeSymmetric(A);
    EXPECT_EQ(IsSymmetric(A),true);

}


TYPED_TEST_P(MatrixTests,OverloadedOperators1)
{
    typedef DMatrix<TypeParam> MatrixT;
    index_t M=10,N=5;
    MatrixT A(M,N),B(A);
    Fill<TypeParam>(A,1.0);
    B=A;
    EXPECT_EQ(A,B);
    EXPECT_EQ(A==B,true);
    EXPECT_EQ(A!=B,false);
    EXPECT_EQ((A+B)(1,1),2.0); // Use Fill constructors EXPECT_EQ((A+B),MatrixT(A.GetLimits(),2.0)
    //EXPECT_EQ((A+B),MatrixT(A.GetLimits(),2.0); TODO
    EXPECT_EQ((A-B   )(1,1),0.0);
    EXPECT_EQ((A*3.0 )(1,1),3.0);
    EXPECT_EQ((A*3   )(1,1),3.0);
    EXPECT_EQ((3.0*A )(1,1),3.0);
    EXPECT_EQ((3*A   )(1,1),3.0);
    EXPECT_EQ((A/2.0 )(1,1),0.5);
    EXPECT_EQ((A/2   )(1,1),0.5);
    EXPECT_EQ((A+=B  )(1,1),2.0);
    EXPECT_EQ((A-=B  )(1,1),1.0);
    EXPECT_EQ((A*=3.0)(1,1),3.0);
    EXPECT_EQ((A*=3  )(1,1),9.0);
    EXPECT_EQ((A/=3.0)(1,1),3.0);
    EXPECT_EQ((A/=3  )(1,1),1.0);

}

TYPED_TEST_P(MatrixTests,MatrixAlgebra)
{
    typedef DMatrix<TypeParam> MatrixT;
    typedef Vector <TypeParam> VectorT;
    index_t M=2,N=4;
    MatrixT A(M,N),B;
    FillLinear<TypeParam>(A,1.0/4,8.0/4);
    B=Transpose(A);
    EXPECT_EQ((A*B).GetLimits(),MatLimits(2,2));
    EXPECT_EQ(Sum( A *B), 25.25);
    EXPECT_EQ(Sum(-A *B),-25.25);
    EXPECT_EQ(Sum( A*-B),-25.25);
    EXPECT_EQ(Sum(-A*-B), 25.25);
    EXPECT_EQ(Sum(A*B-B*A),22.);

    VectorT VL(M),VR(N);
    Fill<TypeParam>(VL, 0.5);
    Fill<TypeParam>(VR,-2.0);
    EXPECT_EQ(VL*A*VR,-9.);
    EXPECT_EQ(VR*B*VL,-9.);
    EXPECT_EQ(Sum(OuterProduct(VL,VR)),-8.);
    EXPECT_EQ(Sum(OuterProduct(VL,VR)*OuterProduct(VR,VL)),16.);
}

TYPED_TEST_P(MatrixTests,MatrixAlgebraPerformance)
{
    typedef DMatrix<TypeParam> MatrixT;
    typedef Vector <TypeParam> VectorT;
#ifdef DEBUG
    index_t N=20,Nreplicates=10;
#else
    int N=200,Nreplicates=10;
#endif
    MatrixT A(N,N),B(N,N),C(N,N);
    VectorT V(N);
    FillRandom<TypeParam>(A);
    FillRandom<TypeParam>(B);
    FillRandom<TypeParam>(V);

    StopWatch s;
    s.Start();
    for (int i=1;i<=Nreplicates;i++) C=A*B;
    s.Stop();
    double Flops=Nreplicates*N*N*N*2/s.GetTime();
    std::cout << "Exp template " << Flops*1e-6 << " MFlops" << std::endl;
    s.Start();
    for (int i=1;i<=Nreplicates;i++) mmul(C,A,B);
    s.Stop();
    Flops=Nreplicates*N*N*N*2/s.GetTime();
    std::cout << "Hand coaded  " << Flops*1e-6 << " MFlops" << std::endl;

    N=1000;
    A.SetLimits(N,N);
    V.SetLimits(N);
    FillRandom<TypeParam>(A);
    FillRandom<TypeParam>(V);
    TypeParam d1, d2=0;
    s.Start();
    for (int i=1;i<=Nreplicates;i++) d1=V*A*V;
    s.Stop();
    double MFlops=(Nreplicates*1e-6)*2*(N*N+N)/s.GetTime();
    std::cout << "Exp template " << MFlops << " MFlops" << std::endl;


    s.Start();
    for (int i=1;i<=Nreplicates;i++) d2=vmul(V,A,V);
    s.Stop();
    assert(fabs(d1-d2)<1e-6);
    MFlops=(Nreplicates*1e-6)*2*(N*N+N)/s.GetTime();
    std::cout << "Hand coaded 1 " << MFlops << " MFlops" << std::endl;

    s.Start();
    for (int i=1;i<=Nreplicates;i++) d2=vmul1(V,A,V);
    s.Stop();
    assert(fabs(d1-d2)<1e-6);
    MFlops=(Nreplicates*1e-6)*2*(N*N+N)/s.GetTime();
    std::cout << "Hand coaded 2 " << MFlops << " MFlops" << std::endl;

}

//inline double fabs(std::complex<double>& c) {return std::abs(c);}

TYPED_TEST_P(MatrixTests,AsciiAndBinaryIO)
{
    typedef DMatrix<TypeParam> MatrixT;
    StreamableObject::SetToAscii();
    {
        MatrixT A(-10,20,-5,50),B;
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
        MatrixT A(-10,20,-5,5),B;
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



template <class T> void mmul(DMatrix<T>& C,const DMatrix<T>& A, const DMatrix<T>& B)
{
  int N=A.GetLimits().Row.High;
  typename DMatrix<T>::Subscriptor s(C);
  for (int i=1;i<=N;i++)
    for (int j=1;j<=N;j++)
    {
      T t(0);
      for (int k=1;k<=N;k++) t+=A(i,k)*B(k,j);
      s(i,j)=t;
    }
}

template <class T> T vmul(const Vector<T> V1,const DMatrix<T>& A, const Vector<T>& V2)
{
  T ret=0;
  int N=A.GetLimits().Row.High;
  for (int i=1;i<=N;i++)
  {
    T t(0);
    for (int j=1;j<=N;j++)
        t+=A(i,j)*V2(j);
    ret+=V1(i)*t;
  }
  return ret;
}
template <class T> T vmul1(const Vector<T> V1,const DMatrix<T>& A, const Vector<T>& V2)
{
  T ret=0;
  int N=A.GetLimits().Row.High;
  for (int j=1;j<=N;j++)
  {
    T t(0);
    for (int i=1;i<=N;i++)
        t+=V1(i)*A(i,j);
    ret+=V2(j)*t;
  }
  return ret;
}

TEST_F(MatrixRealTests,MinMax)
{
    typedef DMatrix<double> MatrixT;
    int M=2,N=4,mn=M*N;
    MatrixT A(M,N),B(A);
    FillLinear(A,1.0/4,2.0); //all elements should be exactly represented in floating point

    EXPECT_EQ(Max(A),2.0);
    EXPECT_EQ(Min(A),0.25);
    {
        MatrixT Ran(100,10);
        FillRandom(Ran,100.0);
        double min=Ran(1,1), max=Ran(1,1);
        for (int i=1; i<=100; i++)
            for (int j=1; j<=10; j++)
            {
                min=Min(min,Ran(i,j));
                max=Max(max,Ran(i,j));
            }
        EXPECT_EQ(Min(Ran),min);
        EXPECT_EQ(Max(Ran),max);
    }

    MatrixT A1(10,5);
    FillRandom(A1);
    MatrixT A2=Transpose(A1);
    EXPECT_EQ(~A2,A1); //Transpose operator
}

TEST_F(MatrixComplexTests,fabsHermitianConj)
{
    typedef std::complex<double> dcmplx;
    typedef DMatrix<dcmplx> MatrixCT;
    typedef Vector <dcmplx> VectorCT;
    typedef DMatrix<double> MatrixRT;
    typedef Vector <double> VectorRT;

    int N=10;
    MatrixCT A(N,N),B(A);
    Fill(A,dcmplx(2,-2));
    B=conj(A);
    EXPECT_EQ(real(Sum(A+B)),400.);
    EXPECT_EQ(Sum(real(A+B)),400.);
    EXPECT_EQ(imag(Sum(A+B)),  0.);
    EXPECT_EQ(Sum(imag(A+B)),  0.);
    EXPECT_EQ(real(Sum(A*B)),8000.);
    EXPECT_EQ(Sum(real(A*B)),8000.);
    EXPECT_EQ(imag(Sum(A*B)),  0.);
    EXPECT_EQ(Sum(imag(A*B)),  0.);

    FillRandom(A);
    B=~A;
    EXPECT_EQ(B==Transpose(conj(A)),true);
    EXPECT_EQ(B(1,2),conj(A(2,1)));
    B=(A+~A)/2.;
    EXPECT_EQ(IsHermitian(B),true);
    VectorCT VL(N),VR(N);
    FillRandom(VL);
    VR=conj(VL);
    dcmplx s=VL*B*VR;
    EXPECT_NEAR(imag(s),0.0,1e-14);
}

TEST_F(MatrixComplexTests,MixedTypes)
{
    typedef std::complex<double> dcmplx;
    typedef DMatrix<dcmplx> MatrixCT;
    typedef Vector <dcmplx> VectorCT;
    typedef DMatrix<double> MatrixRT;
    typedef Vector <double> VectorRT;

    int N=10;
    MatrixCT Ac(N,N);
    MatrixRT Ar(N,N);
    Fill(Ac,dcmplx(2,-2));
    Fill(Ar,0.5);
    VectorCT Vc(N,dcmplx(0.25,0.5));
    VectorRT Vr(N,2.0);
    dcmplx   Sc(0.5,-0.25);
    double   Sr(4.0);
    EXPECT_EQ((Ac*Ar)(1,2),dcmplx(10,-10));
    EXPECT_EQ((Ar*Ac)(1,2),dcmplx(10,-10));
    EXPECT_EQ((Ac*Vr)(1),dcmplx(40,-40));
    EXPECT_EQ((Vr*Ac)(1),dcmplx(40,-40));
    EXPECT_EQ((Vc*Ar)(1),dcmplx(1.25,2.5));
    EXPECT_EQ((Ar*Vc)(1),dcmplx(1.25,2.5));
    EXPECT_EQ((Vc*Vc)     ,dcmplx(-1.875,2.5));
    EXPECT_EQ((Vc*Vr)     ,dcmplx(5,10));
    EXPECT_EQ((Vr*Vc)     ,dcmplx(5,10));
    EXPECT_EQ((Ac*Sr)(1,2),dcmplx(8,-8));
    EXPECT_EQ((Sr*Ac)(1,2),dcmplx(8,-8));
    EXPECT_EQ((Sc*Ar)(1,2),dcmplx(.25,-.125));
    EXPECT_EQ((Ar*Sc)(1,2),dcmplx(.25,-.125));
    EXPECT_EQ((Vc*Sr)(1),dcmplx(1,2));
    EXPECT_EQ((Sr*Vc)(1),dcmplx(1,2));
    EXPECT_EQ((Sc*Vr)(1),dcmplx(1,-0.5));
    EXPECT_EQ((Vr*Sc)(1),dcmplx(1,-0.5));
}


REGISTER_TYPED_TEST_SUITE_P(
            MatrixTests,
            Constructors,
            Fill_SetLimits_SubMatrix,
            SumDotMaxMinDirectMultiply,
            Transpose_Slices,
            OverloadedOperators1,
            MatrixAlgebra,
            MatrixAlgebraPerformance,
            AsciiAndBinaryIO
            );
using MyTypes = ::testing::Types<double,std::complex<double>>;
INSTANTIATE_TYPED_TEST_SUITE_P(My, MatrixTests, MyTypes);
