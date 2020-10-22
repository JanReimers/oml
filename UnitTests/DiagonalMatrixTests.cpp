#include "oml/diagonalmatrix.h"
#include "oml/matrix.h"
#include "oml/vector.h"
#include "oml/random.h"
#include "stopw.h"
#include "gtest/gtest.h"
#include <iostream>
#include <fstream>
#include <complex>

using std::cout;
using std::endl;

template <class T> class DiagonalMatrixTests : public ::testing::Test
{
public:
    DiagonalMatrixTests()
    {
        StreamableObject::SetToPretty();
    }
};

class DiagonalMatrixComplexTests : public ::testing::Test
{
public:
    DiagonalMatrixComplexTests()
    {
        StreamableObject::SetToPretty();
    }
};

TYPED_TEST_SUITE_P(DiagonalMatrixTests);

TYPED_TEST_P(DiagonalMatrixTests,Constructors)
{
    typedef DiagonalMatrix<TypeParam> MatrixT;
    typedef Vector<TypeParam> VectorT;

    MatrixT A0;
    EXPECT_EQ(A0.GetNumRows(),0);
    EXPECT_EQ(A0.GetNumCols(),0);

    MatrixT A1(10);
    EXPECT_EQ(A1.GetNumRows(),10);
    EXPECT_EQ(A1.GetNumCols(),10);

    MatLimits lim(0,5,0,5);
    MatrixT A2(lim);
    EXPECT_EQ(A2.GetNumRows(),6);
    EXPECT_EQ(A2.GetNumCols(),6);
    EXPECT_EQ(A2.GetLimits(),lim);
    EXPECT_EQ(A2.GetRowLow(),0);
    EXPECT_EQ(A2.GetColLow(),0);
    EXPECT_EQ(A2.size(),6);

    MatrixT A3(A2),A4,A5(4);
    EXPECT_EQ(A3.GetLimits(),A2.GetLimits());
    EXPECT_NE(A4.GetLimits(),A2.GetLimits());
    EXPECT_NE(A5.GetLimits(),A2.GetLimits());
    A4=A2;
    A5=A2;
    EXPECT_EQ(A4.GetLimits(),A2.GetLimits());
    EXPECT_EQ(A5.GetLimits(),A2.GetLimits());

    VectorT V(7);
    FillRandom(V);
    MatrixT A6(V);
    A1=V;
    EXPECT_TRUE(A1==A6);
    EXPECT_EQ(A1.GetDiagonal(),V);

}

TYPED_TEST_P(DiagonalMatrixTests,Fill_SetLimits_SubMatrix)
{
    typedef DiagonalMatrix<TypeParam> MatrixT;
    MatrixT A1(5);
    FillRandom(A1);
    MatLimits newLimits(2,4,2,4);
    MatrixT A2=A1.SubMatrix(newLimits);
    A1.SetLimits(newLimits,true);
    EXPECT_EQ(A1,A2);
}

TYPED_TEST_P(DiagonalMatrixTests,Slices)
{
    typedef DiagonalMatrix<TypeParam> MatrixT;
    typedef Vector <TypeParam> VectorT;

    index_t N=10;
    MatrixT A3(N);
    FillRandom(A3);
    VectorT D=A3.GetDiagonal();
    EXPECT_EQ(D.size(),N);
    for (index_t i=1; i<=N; i++)
        EXPECT_EQ(A3(i,i),D(i));
}

TYPED_TEST_P(DiagonalMatrixTests,SumDotMaxMinDirectMultiply)
{
    typedef DiagonalMatrix<TypeParam> MatrixT;
    index_t N=8;
    MatrixT A(N),B(A);
    FillLinear<TypeParam>(A,1.0/4,2.0); //all elements should be exactly represented in floating point
    EXPECT_EQ(Sum(A),N*(N+1)/2/4.);
    EXPECT_EQ(Dot(A,A),51.0/4); //Look up sum n^2 rule
    Fill<TypeParam>(B,0.5);
    EXPECT_EQ(DirectMultiply(A,B)(1,1),1.0/8);
    EXPECT_EQ(DirectDivide(A,B)(1,1),1.0/2);

    A.GetDiagonal().Fill(FillType::Unit);
    EXPECT_EQ(A(2,2),1.0);
    EXPECT_EQ(A(2,1),0.0);
    EXPECT_EQ(Sum(A),TypeParam(N));

}

TYPED_TEST_P(DiagonalMatrixTests,OverloadedOperators1)
{
    typedef DiagonalMatrix<TypeParam> MatrixT;
    index_t N=10;
    MatrixT A(N),B(A);
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

TYPED_TEST_P(DiagonalMatrixTests,MatrixAlgebra)
{
    typedef DiagonalMatrix<TypeParam> MatrixT;
    typedef Vector <TypeParam> VectorT;
    index_t N=8;
    MatrixT A(N),B;
    FillLinear<TypeParam>(A,1.0/4,8.0/4);
    TypeParam sd=Sum(A.GetDiagonal());
    B=A*0.5;
    EXPECT_EQ((A*B).GetLimits(),MatLimits(N,N));
    EXPECT_EQ(Sum( A *B), 6.375);
    EXPECT_EQ(Sum(-A *B),-6.375);
    EXPECT_EQ(Sum( A*-B),-6.375);
    EXPECT_EQ(Sum(-A*-B), 6.375);
    EXPECT_EQ(Sum(A*B-B*A),0.0);

    VectorT VL(N),VR(N);
    Fill<TypeParam>(VL, 0.5);
    Fill<TypeParam>(VR,-2.0);
    EXPECT_EQ(VL*A*VR,-sd);
    EXPECT_EQ(VR*B*VL,-sd/2.0);
    B=A*16;
    MatrixT C=A*B*B*A*A*A*A*B;
    EXPECT_EQ(C(1,1), 0.0625);

    Matrix<TypeParam> F(A.GetLimits());
    Fill(F,TypeParam(2));
    EXPECT_EQ((A*F)(1,1), 0.5);
    EXPECT_EQ((F*A)(1,1), 0.5);
    EXPECT_EQ((F*A*A*F*F*F*A*F*A)(1,1), 58752.);

}

TEST_F(DiagonalMatrixComplexTests,fabsHermitianConj)
{
    typedef std::complex<double> dcmplx;
    typedef DiagonalMatrix<dcmplx> MatrixCT;
    typedef Vector        <dcmplx> VectorCT;
    typedef DiagonalMatrix<double> MatrixRT;
    typedef Vector        <double> VectorRT;

    int N=10;
    MatrixCT A(N),B(A);
    Fill(A,dcmplx(2,-2));
    B=conj(A);
    EXPECT_EQ(real(Sum(A+B)),40.);
    EXPECT_EQ(Sum(real(A+B)),40.);
    EXPECT_EQ(imag(Sum(A+B)), 0.);
    EXPECT_EQ(Sum(imag(A+B)), 0.);
    EXPECT_EQ(real(Sum(A*B)),80.);
    EXPECT_EQ(Sum(real(A*B)),80.);
    EXPECT_EQ(imag(Sum(A*B)), 0.);
    EXPECT_EQ(Sum(imag(A*B)), 0.);

}

TEST_F(DiagonalMatrixComplexTests,MixedTypes)
{
    typedef std::complex<double> dcmplx;
    typedef DiagonalMatrix<dcmplx> MatrixCT;
    typedef        Vector <dcmplx> VectorCT;
    typedef DiagonalMatrix<double> MatrixRT;
    typedef        Vector <double> VectorRT;

    int N=10;
    MatrixCT Ac(N);
    MatrixRT Ar(N);
    Fill(Ac,dcmplx(2,-2));
    Fill(Ar,0.5);
    VectorCT Vc(N,dcmplx(0.25,0.5));
    VectorRT Vr(N,2.0);
    dcmplx   Sc(0.5,-0.25);
    double   Sr(4.0);
    // Try every combination we can think of
    EXPECT_EQ((Ac*Ar)(2,2),dcmplx(1,-1));
    EXPECT_EQ((Ar*Ac)(2,2),dcmplx(1,-1));
    EXPECT_EQ((Ac*Vr)(1)  ,dcmplx(4,-4));
    EXPECT_EQ((Vr*Ac)(1)  ,dcmplx(4,-4));
    EXPECT_EQ((Vc*Ar)(1)  ,dcmplx(.125,.25));
    EXPECT_EQ((Ar*Vc)(1)  ,dcmplx(.125,.25));
    EXPECT_EQ((Vc*Vc)     ,dcmplx(-1.875,2.5));
    EXPECT_EQ((Vc*Vr)     ,dcmplx(5,10));
    EXPECT_EQ((Vr*Vc)     ,dcmplx(5,10));
    EXPECT_EQ((Ac*Sr)(2,2),dcmplx(8,-8));
    EXPECT_EQ((Sr*Ac)(2,2),dcmplx(8,-8));
    EXPECT_EQ((Sc*Ar)(2,2),dcmplx(.25,-.125));
    EXPECT_EQ((Ar*Sc)(2,2),dcmplx(.25,-.125));
    EXPECT_EQ((Vc*Sr)(1)  ,dcmplx(1,2));
    EXPECT_EQ((Sr*Vc)(1)  ,dcmplx(1,2));
    EXPECT_EQ((Sc*Vr)(1)  ,dcmplx(1,-0.5));
    EXPECT_EQ((Vr*Sc)(1)  ,dcmplx(1,-0.5));
    // Try some long chains of mixed type multiplications
    EXPECT_EQ((Ac*Ar*Ac*Ar*Ar*Ac)(2,2),dcmplx(-2,-2));
    EXPECT_EQ((Vc*Ac*Ar*Ac*Ar*Ar*Ac*Vr),dcmplx(10,-30));
    EXPECT_EQ((Ac*Ar*Sc*Ar*Ar*Sc)(2,2),dcmplx(-.015625,-.109375));
    EXPECT_EQ((Vc*Ac*Sc*Ar*Sc*Ar*Ar*Ac*Vr),dcmplx(.62500,-3.43750));

    EXPECT_EQ((Ac+Sc)(1,1),dcmplx( 2.5,-2.25));
    EXPECT_EQ((Ac-Sc)(1,1),dcmplx( 1.5,-1.75));
    EXPECT_EQ((Sc+Ac)(1,1),dcmplx( 2.5,-2.25));
    EXPECT_EQ((Sc-Ac)(1,1),dcmplx(-1.5, 1.75));
    EXPECT_EQ((Ac+Sr)(1,1),dcmplx( 6.0,-2.0 ));
    EXPECT_EQ((Ac-Sr)(1,1),dcmplx(-2.0,-2.0 ));
    EXPECT_EQ((Sr+Ac)(1,1),dcmplx( 6.0,-2.0 ));
    EXPECT_EQ((Sr-Ac)(1,1),dcmplx( 2.0, 2.0 ));
    EXPECT_EQ((Ar+Sr)(1,1), 4.5);
    EXPECT_EQ((Ar-Sr)(1,1),-3.5);
    EXPECT_EQ((Sr+Ar)(1,1), 4.5);
    EXPECT_EQ((Sr-Ar)(1,1), 3.5);
    EXPECT_EQ((Ar+Sc)(1,1),dcmplx( 1.0,-.25 ));
    EXPECT_EQ((Ar-Sc)(1,1),dcmplx( 0.0, .25 ));
    EXPECT_EQ((Sc+Ar)(1,1),dcmplx( 1.0,-.25 ));
    EXPECT_EQ((Sc-Ar)(1,1),dcmplx( 0.0,-.25 ));
    EXPECT_EQ((Ar+1)(1,1), 1.5);
    EXPECT_EQ((Ar-1)(1,1),-0.5);
    EXPECT_EQ((1+Ar)(1,1), 1.5);
    EXPECT_EQ((1-Ar)(1,1), 0.5);

    N=2;
    Matrix<dcmplx> U(N,N),Vdagger(N,N),Ac1(N,N);
    Ar.SetLimits(N);
    Fill(U,dcmplx(1.0));
    Fill(Vdagger,dcmplx(1.0));
//    Fill(Ar,0.5);
    Ar(1)=2.0;
    Ar(2)=0.5;

    Fill(Ac1,dcmplx(1.0));
    Matrix<dcmplx> R=U*Ar*Vdagger; //Make sure this compiles
//    Matrix<dcmplx> R=U*Ar; //Make sure this compiles
    cout << R << endl;
}




REGISTER_TYPED_TEST_SUITE_P(DiagonalMatrixTests,
                            Constructors,
                            Fill_SetLimits_SubMatrix,
                            Slices,
                            SumDotMaxMinDirectMultiply,
                            OverloadedOperators1,
                            MatrixAlgebra
                            );
//using MyTypes = ::testing::Types<double>;
using MyTypes = ::testing::Types<double,std::complex<double>>;
INSTANTIATE_TYPED_TEST_SUITE_P(My, DiagonalMatrixTests, MyTypes);

