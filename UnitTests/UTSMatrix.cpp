// File: UTSMatrix.cpp  Unit test the SMatrix<T> class for double/complex data types.

//#include "stopw.h"
#include "gtest/gtest.h"
// #include <omp.h>
#include <iostream>
#include <iomanip>
#include <fstream>
#include <complex>

import oml.Vector;
import oml.Matrix;
import oml.SMatrix;
import omlNumericalRecipes;



template <class T> class SMatrixTests : public ::testing::Test
{
public:
    SMatrixTests()
    {
        StreamableObject::SetToPretty();
    }
};

class SMatrixDoubleTests : public ::testing::Test
{
public:
    SMatrixDoubleTests()
    {
        StreamableObject::SetToPretty();
    }
};

class SMatrixComplexTests : public ::testing::Test
{
public:
    SMatrixComplexTests()
    {
        StreamableObject::SetToPretty();
    }
};

TYPED_TEST_SUITE_P(SMatrixTests);

TYPED_TEST_P(SMatrixTests,Constructors)
{
    typedef SMatrix<TypeParam> SMat;
    typedef Matrix<TypeParam> Mat;

    SMat A0;
    EXPECT_EQ(A0.GetNumRows(),0);
    EXPECT_EQ(A0.GetNumCols(),0);

    SMat A1(10);
    EXPECT_EQ(A1.GetNumRows(),10);
    EXPECT_EQ(A1.GetNumCols(),10);

    MatLimits lim(0,5,0,5);
    SMat A2(lim);
    EXPECT_EQ(A2.GetNumRows(),6);
    EXPECT_EQ(A2.GetNumCols(),6);
    EXPECT_EQ(A2.GetLimits(),lim);
    EXPECT_EQ(A2.GetRowLow(),0);
    EXPECT_EQ(A2.GetColLow(),0);
    EXPECT_EQ(A2.size(),6*7/2);

    SMat A3(A2),A4,A5(4);
    EXPECT_EQ(A3.GetLimits(),A2.GetLimits());
    EXPECT_NE(A4.GetLimits(),A2.GetLimits());
    EXPECT_NE(A5.GetLimits(),A2.GetLimits());
    A4=A2;
    A5=A2;
    EXPECT_EQ(A4.GetLimits(),A2.GetLimits());
    EXPECT_EQ(A5.GetLimits(),A2.GetLimits());
    //
    //  Cross contruction
    //
    FillRandom(A1);
    A1.GetDiagonal()=A1.GetDiagonal()+conj(A1.GetDiagonal()); //Hermitian.
    Mat B(A1); //Easy
    SMat C(B); //B check for symmetry.
    EXPECT_EQ(A1,C);
}

TYPED_TEST_P(SMatrixTests,Fill_SetLimits_SubMatrix)
{
    typedef SMatrix<TypeParam> MatrixT;
    MatrixT A1(5);
    FillRandom(A1);
    MatLimits newLimits(2,4,2,4);
    MatrixT A2=A1.SubMatrix(newLimits);
    A1.SetLimits(newLimits,true);
    EXPECT_EQ(A1,A2);
}
TYPED_TEST_P(SMatrixTests,Slices)
{
    typedef SMatrix<TypeParam> MatrixT;
    typedef Vector <TypeParam> VectorT;

    index_t N=10;
    MatrixT A3(N);
    FillRandom(A3);
    VectorT D=A3.GetDiagonal();
    EXPECT_EQ(D.size(),N);
    for (index_t i=1; i<=N; i++)
        EXPECT_EQ(A3(i,i),D(i));
}

TYPED_TEST_P(SMatrixTests,SumDotMaxMinDirectMultiply)
{
    typedef SMatrix<TypeParam> MatrixT;
    typedef  Matrix<TypeParam> MatrixFT;
    index_t N=4;
    MatrixT A(N),B(A);
    FillLinear<TypeParam>(A,1.0/8,10.0/8); //all elements should be exactly represented in floating point
    MatrixFT Af(A),Bf(B);
    EXPECT_EQ(Sum(A),Sum(Af));
    EXPECT_EQ(Dot(A,A),Dot(Af,Af)); //Look up sum n^2 rule
    Fill<TypeParam>(B,0.5);
    Bf=B;
    EXPECT_EQ(DirectMultiply(A,B),DirectMultiply(Af,Bf));
    EXPECT_EQ(DirectDivide  (A,B),DirectDivide  (Af,Bf));

    Unit(A);
    TypeParam a12=A(1,2);
    EXPECT_EQ(A(1,2),a12);
    EXPECT_EQ(A(2,2),1.0);
    EXPECT_EQ(A(1,2),0.0);
    EXPECT_EQ(Sum(A),TypeParam(N));

}

TYPED_TEST_P(SMatrixTests,OverloadedOperators1)
{
    typedef SMatrix<TypeParam> MatrixT;
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

TYPED_TEST_P(SMatrixTests,MatrixAlgebra)
{
    typedef Matrix<TypeParam>  MatrixFT;
    typedef SMatrix<TypeParam> MatrixT;
    typedef Vector <TypeParam> VectorT;
    index_t N=4;
    MatrixT A(N),B;
    FillLinear<TypeParam>(A,1.0/8,10.0/8);
    TypeParam sd=Sum(A);
    B=A*0.5;
    MatrixFT Af(A),Bf(B);
    TypeParam s=Sum( Af *Bf);
    EXPECT_EQ((A*B).GetLimits(),MatLimits(N,N));
    EXPECT_EQ(Sum( A *B), s);
    EXPECT_EQ(Sum(-A *B),-s);
    EXPECT_EQ(Sum( A*-B),-s);
    EXPECT_EQ(Sum(-A*-B), s);
    EXPECT_EQ(Sum(A*B-B*A),0.0);

    VectorT VL(N),VR(N);
    Fill<TypeParam>(VL, 0.5);
    Fill<TypeParam>(VR,-2.0);
    EXPECT_EQ(VL*A*VR,-sd);
    EXPECT_EQ(VR*B*VL,-sd/2.0);
    B=A*2;
    MatrixFT C=A*B*B*A*A*A*A*B;
    Af=A;Bf=B;
    MatrixFT Cf=Af*Bf*Bf*Af*Af*Af*Af*Bf;
    EXPECT_EQ(C,Cf);

    Matrix<TypeParam> F(A.GetLimits()),F2;
    Fill(F,TypeParam(2));
    EXPECT_EQ((A*F)(1,1), 2.5);
    EXPECT_EQ((F*A)(1,1), 2.5);
    EXPECT_EQ((F*A*A*F*F*F*A*F*A)(1,1), 223170.);
    EXPECT_EQ((A+F)(1,1), 2.125);
    EXPECT_EQ((F+A)(1,1), 2.125);
    EXPECT_EQ((A-F)(1,1),-1.875);
    EXPECT_EQ((F-A)(1,1), 1.875);
    F2=F+A;
    //B=F+A; This should not compile Symmetric = Full + Symmetric, no way store the result in Symmetric.
    F2=F*A;
    F2*=A;
}


TEST_F(SMatrixComplexTests,fabsHermitianConj)
{
    typedef std::complex<double> dcmplx;
    typedef SMatrix<dcmplx> MatrixCT;
    typedef  Matrix<dcmplx> MatrixFCT;

    int N=3;
    MatrixCT A(N),B(A);
    MatrixFCT Af(N,N),Bf(Af);
    Fill(A,dcmplx(2,-2));
    Af=A;
    B=conj(A);
    Bf=B;
    EXPECT_EQ(A+B,Af+Bf);
    EXPECT_EQ(A-B,Af-Bf);
    EXPECT_EQ(A*B,Af*Bf);
    EXPECT_EQ(real(Sum(A+B)),real(Sum(Af+Bf)));
    EXPECT_EQ(Sum(real(A+B)),Sum(real(Af+Bf)));
    EXPECT_EQ(imag(Sum(A+B)),imag(Sum(Af+Bf)));

    EXPECT_EQ(Sum(imag(A+B)),Sum(imag(Af+Bf)));
    EXPECT_EQ(real(Sum(A*B)),real(Sum(Af*Bf)));
    EXPECT_EQ(Sum(real(A*B)),Sum(real(Af*Bf)));
    EXPECT_EQ(imag(Sum(A*B)),imag(Sum(Af*Bf)));
    EXPECT_EQ(Sum(imag(A*B)),Sum(imag(Af*Bf)));

}

void FillPos(SMatrix<double>& m)
{
    typename SMatrix<double>::Subscriptor s(m);
    for (int i=m.GetLimits().Row.Low; i<=m.GetLimits().Row.High; i++)
        for (int j=m.GetLimits().Col.Low; j<=m.GetLimits().Col.High; j++)
        {
            double del=(i-j)/2.0;
            s(i,j)=exp(-del*del);
        }
}

TEST_F(SMatrixDoubleTests,LinearAlgebra)
{
    typedef Matrix<double>  MatrixFT;
    typedef SMatrix<double> MatrixT;
    typedef Vector <double> VectorT;
    index_t N=10;
    MatrixT A(N);
    VectorT ones(N);
    FillPos(A);
    Fill(ones,1.0);
        
    MatrixT Ai=InvertSymmetric(A);
    MatrixFT I1=A*Ai,I2=Ai*A;
    I1.GetDiagonal()=I1.GetDiagonal()-ones;
    I2.GetDiagonal()=I2.GetDiagonal()-ones;
    EXPECT_NEAR(Max(fabs(I1)),0,1e-12);
    EXPECT_NEAR(Max(fabs(I2)),0,1e-12);
}

TEST_F(SMatrixComplexTests,MixedTypes)
{
    typedef std::complex<double> dcmplx;
    typedef SMatrix<dcmplx> MatrixCT;
    typedef  Matrix<dcmplx> MatrixFCT;
    typedef  Vector<dcmplx> VectorCT;
    typedef SMatrix<double> MatrixRT;
    typedef  Matrix<double> MatrixFRT;
    typedef  Vector<double> VectorRT;

    int N=10;
    MatrixCT Ac(N);
    MatrixRT Ar(N);
    Fill(Ac,dcmplx(2,-2));
    Fill(Ar,0.5);
    MatrixFCT Acf(Ac);
    MatrixFRT Arf(Ar);



    VectorCT Vc(N,dcmplx(0.25,0.5));
    VectorRT Vr(N,2.0);
    dcmplx   Sc(0.5,-0.25);
    double   Sr(4.0);
    // Try every combination we can think of
    EXPECT_EQ(Ac*Ar,Acf*Arf);
    EXPECT_EQ(Ar*Ac,Arf*Acf);
    EXPECT_EQ(Ac*Vr,Acf*Vr );
    EXPECT_EQ(Vr*Ac,Vr *Acf);
    EXPECT_EQ(Vc*Ar,Vc *Arf);
    EXPECT_EQ(Ar*Vc,Arf*Vc );
    EXPECT_EQ(Ac*Sr,Acf*Sr );
    EXPECT_EQ(Sr*Ac,Sr *Acf);
    EXPECT_EQ(Sc*Ar,Sc *Arf);
    EXPECT_EQ(Ar*Sc,Arf*Sc );
    // Try some long chains of mixed type multiplications
    EXPECT_EQ(Ac*Ar*Ac*Ar*Ar*Ac,Acf*Arf*Acf*Arf*Arf*Acf);
    EXPECT_EQ(Vc*Ac*Ar*Ac*Ar*Ar*Ac*Vr,Vc*Acf*Arf*Acf*Arf*Arf*Acf*Vr);
    EXPECT_EQ(Ac*Ar*Sc*Ar*Ar*Sc,Acf*Arf*Sc*Arf*Arf*Sc);
    EXPECT_EQ(Vc*Ac*Sc*Ar*Sc*Ar*Ar*Ac*Vr,Vc*Acf*Sc*Arf*Sc*Arf*Arf*Acf*Vr);

    EXPECT_EQ(Ac+Sc,Acf+Sc );
    EXPECT_EQ(Ac-Sc,Acf-Sc );
    EXPECT_EQ(Sc+Ac,Sc +Acf);
    EXPECT_EQ(Sc-Ac,Sc -Acf);
    EXPECT_EQ(Ac+Sr,Acf+Sr );
    EXPECT_EQ(Ac-Sr,Acf-Sr );
    EXPECT_EQ(Sr+Ac,Sr +Acf);
    EXPECT_EQ(Sr-Ac,Sr -Acf);
    EXPECT_EQ(Ar+Sr,Arf+Sr );
    EXPECT_EQ(Ar-Sr,Arf-Sr );
    EXPECT_EQ(Sr+Ar,Sr +Arf);
    EXPECT_EQ(Sr-Ar,Sr -Arf);
    EXPECT_EQ(Ar+Sc,Arf+Sc );
    EXPECT_EQ(Ar-Sc,Arf-Sc );
    EXPECT_EQ(Sc+Ar,Sc +Arf);
    EXPECT_EQ(Sc-Ar,Sc -Arf);
    EXPECT_EQ(Ar+1, Arf+1);
    EXPECT_EQ(Ar-1, Arf-1);
    EXPECT_EQ(1+Ar, 1+Arf);
    EXPECT_EQ(1-Ar, 1-Arf);

}



REGISTER_TYPED_TEST_SUITE_P(SMatrixTests
                            ,Constructors
                            ,Fill_SetLimits_SubMatrix
                            ,Slices
                            ,SumDotMaxMinDirectMultiply
                            ,OverloadedOperators1
                            ,MatrixAlgebra
                            );
// using MyTypes = ::testing::Types<double>;
using MyTypes = ::testing::Types<double,std::complex<double>>;
INSTANTIATE_TYPED_TEST_SUITE_P(My, SMatrixTests, MyTypes);



