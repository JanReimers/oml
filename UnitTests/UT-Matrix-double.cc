// File: UT-Matrix-double.cc  Unit test the Matrix class for double data types.

// Copyright (1994-2005), Jan N. Reimers

#include "oml/UnitTests/UnitTest.H"
#include "oml/vector_io.h"
#include "oml/matrix.h"
#include "oml/smatrix.h"
#include "oml/numeric.h"
#include "oml/lusolver.h"
#include "oml/svsolver.h"
#include "oml/minmax.h"
#include "oml/random.h"
#include "oml/stopw.h"
#include <fstream>
#include <iostream>

using std::cout;
using std::endl;

typedef Matrix <double> MatrixType;
typedef Vector <double> VectorType;


void mmul(MatrixType& C,const MatrixType& A, const MatrixType&B);
double vmul(const VectorType V1,const MatrixType& A, const VectorType& V2);
double vmul1(const VectorType V1,const MatrixType& A, const VectorType& V2);

template <class T, class A, Store M, Data D> void
FillPos(Indexable<T,A,M,D,MatrixShape>& m)
{
    typename A::Subscriptor s(m);
    for (int i=m.GetLimits().Row.Low; i<=m.GetLimits().Row.High; i++)
        for (int j=m.GetLimits().Col.Low; j<=m.GetLimits().Col.High; j++)
        {
            double del=(i-j)/2.0;
            s(i,j)=exp(-del*del);
        }
}

bool TestM1();
bool TestM2();
bool TestM3();
bool TestM4();
bool TestM5();
bool TestM6();
bool TestM7();
bool TestM8();
bool TestM9();

int TestMatrixDouble()
{
    const char* Class="Matrix<double>";
    bool pass=true;
//
//  Start message.
//

    StartClass(Class);
    StreamableObject::SetToPretty();
    std::cout << "Matrix<double> foot print=" << sizeof(Matrix<double>) << " bytes" << std::endl;

    bool p[9]= {true,true,true,true,true,true,true,true,true};

//    p[0]=TestM1();
    p[1]=TestM2();
//    p[2]=TestM3();
//    p[3]=TestM4();
//    p[4]=TestM5();
//    p[5]=TestM6();
//    p[6]=TestM7();
//    p[7]=TestM8();
//    p[8]=TestM9();

    pass=p[0]&&p[1]&&p[2]&&p[3]&&p[4]&&p[5]&&p[6]&&p[7]&&p[8];




    return pass ? 0 : -1;
}



void mmul(MatrixType& C,const MatrixType& A, const MatrixType&B)
{
    int N=A.GetLimits().Row.High;
    MatrixType::Subscriptor s(C);
    for (int i=1; i<=N; i++)
        for (int j=1; j<=N; j++)
        {
            double t=0;
            for (int k=1; k<=N; k++) t+=A(i,k)*B(k,j);
            s(i,j)=t;
        }
}

double vmul1(const VectorType V1,const MatrixType& A, const VectorType& V2)
{
    return V1*A*V1;
}

double vmul(const VectorType V1,const MatrixType& A, const VectorType& V2)
{
    double ret=0;
    int N=A.GetLimits().Row.High;
    for (int i=1; i<=N; i++)
    {
        double t=0;
        for (int j=1; j<=N; j++) t+=A(i,j)*V2(j);
        ret+=V1(i)*t;
    }
    return ret;
}


//###########################################################################
//
//  Constructors
//
bool TestM1()
{
    bool pass=true;

    EXPECT(MatrixType(),"(1:0),(1:0) ");
    EXPECT(MatrixType(2,3),"(1:2),(1:3) \n[ * * * ]\n[ * * * ]\n");
    EXPECT(MatrixType(-2,1,-1,3),
           "(-2:1),(-1:3) \n[ * * * * * ]\n[ * * * * * ]\n[ * * * * * ]\n[ * * * * * ]\n");
    EXPECT(MatrixType(VecLimits(-2,1),VecLimits(-1,3)),
           "(-2:1),(-1:3) \n[ * * * * * ]\n[ * * * * * ]\n[ * * * * * ]\n[ * * * * * ]\n");
    EXPECT(MatrixType(MatLimits(-2,1,-1,3)),
           "(-2:1),(-1:3) \n[ * * * * * ]\n[ * * * * * ]\n[ * * * * * ]\n[ * * * * * ]\n");
    MatrixType A1(3,3),A2(3,3);
    Fill(A1,3.0);
    EXPECT3(temp,MatrixType temp(A1),"(1:3),(1:3) \n[ 3 3 3 ]\n[ 3 3 3 ]\n[ 3 3 3 ]\n");
    EXPECT(A2=A1,"(1:3),(1:3) \n[ 3 3 3 ]\n[ 3 3 3 ]\n[ 3 3 3 ]\n");

    EXPECT2(A1,SetLimits(0,4,0,4,true),
            "(0:4),(0:4) \n[ * * * * * ]\n[ * 3 3 3 * ]\n[ * 3 3 3 * ]\n[ * 3 3 3 * ]\n[ * * * * * ]\n");
    EXPECT2(A1,SetLimits(0,4,0,4),
            "(0:4),(0:4) \n[ * * * * * ]\n[ * * * * * ]\n[ * * * * * ]\n[ * * * * * ]\n[ * * * * * ]\n");

    A1.SetLimits(3,3);
    EXPECT3(A2,FillLinear(A2,1.0,9.0),"(1:3),(1:3) \n[ 1 4 7 ]\n[ 2 5 8 ]\n[ 3 6 9 ]\n");
    EXPECT(A2(2,3),"8");
    A2(3,2)=888;
    EXPECT3(A2,A2(3,2)=888,"(1:3),(1:3) \n[ 1 4 7 ]\n[ 2 5 8 ]\n[ 3 888 9 ]\n");
    std::cout << A1.GetLimits() << "   " << A2.GetLimits() << std::endl;
    EXPECT3(A1,A1=Transpose(A2),"(1:3),(1:3) \n[ 1 2 3 ]\n[ 4 5 888 ]\n[ 7 8 9 ]\n");
    A2.SetLimits(1,3,1,4,true);
    A1.SetLimits(1,4,1,3);
    EXPECT3(A1,A1=Transpose(A2),"(1:4),(1:3) \n[ 1 2 3 ]\n[ 4 5 888 ]\n[ 7 8 9 ]\n[ * * * ]\n");
    A2.SetLimits(3,3,true);
    A1.SetLimits(3,3);
    EXPECT3(A1,A1=Transpose(A2+A2),"(1:3),(1:3) \n[ 2 4 6 ]\n[ 8 10 1776 ]\n[ 14 16 18 ]\n");
    EXPECT3(A1,A1=1.0+Transpose(A2+A2),"(1:3),(1:3) \n[ 3 5 7 ]\n[ 9 11 1777 ]\n[ 15 17 19 ]\n");

    EXPECT(A1.SubMatrix(MatLimits(2,3,1,2)),"(2:3),(1:2) \n[ 9 11 ]\n[ 15 17 ]\n");
    return pass;
}


bool TestM2()
{
    bool pass=true;

    MatrixType A(2,8,2,8);
    FillLinear(A,1.0,49.0);
    EXPECT(A.GetDiagonal(),"(2:8){ 1 9 17 25 33 41 49 }");
    EXPECT(A.GetColumn (3),"(2:8){ 8 9 10 11 12 13 14 }");
    EXPECT(A.GetRow    (6),"(2:8){ 5 12 19 26 33 40 47 }");
    EXPECT(Sum(A.GetColumn (3)),77);
    EXPECT(Sum(A.GetRow    (6)),182);
    EXPECT(Max(A.GetRow    (6)),47);
    EXPECT(Min(A.GetRow    (6)),5);
    EXPECT(Sum(A.GetDiagonal()),175);
    EXPECT(A.GetRow(6)*A.GetColumn (3),2198);
    EXPECT(Min(Vector<double>(A.GetRow    (6))),5.0);
    EXPECT(A.GetRow(7)=A.GetRow(6),"(2:8){ 5 12 19 26 33 40 47 }");
    EXPECT(A.GetRow(8)=Vector<double>(A.GetRow(6)),"(2:8){ 5 12 19 26 33 40 47 }");
    EXPECT(A.GetRow(8)+=A.GetRow(6),"(2:8){ 10 24 38 52 66 80 94 }");
//
// Treat rows and columns as arrays
//
    int colLow=A.GetLimits().Col.Low;
    int rowLow=A.GetLimits().Row.Low;
    Array<double> Arr1=A.GetRowAsArray(2);
    Array<double> Arr2=A.GetColumnAsArray(8);
    for (int i=0; i<Arr1.size(); i++)
    {
        EXPECT(Arr1[i],A(2,i+colLow));
    }
    for (int i=0; i<Arr2.size(); i++)
    {
        EXPECT(Arr2[i],A(i+rowLow,8));
    }


    FillRandom(Arr1);
    FillRandom(Arr2);
    A.GetRowAsArray(2)=Arr1;

    for (int i=0; i<Arr1.size(); i++)
    {
        EXPECT(Arr1[i],A(2,i+colLow));
    }
    A.GetColumnAsArray(8)=Arr2;
    for (int i=0; i<Arr2.size(); i++)
    {
        EXPECT(Arr2[i],A(i+rowLow,8));
    }

    return pass;

}

bool TestM3()
{
    bool pass=true;

    MatrixType A(3,3),B(3,3);
    Fill(A,13.0);
    Fill(B, 5.0);

    EXPECT(A,"(1:3),(1:3) \n[ 13 13 13 ]\n[ 13 13 13 ]\n[ 13 13 13 ]\n");
    EXPECT(B,"(1:3),(1:3) \n[ 5 5 5 ]\n[ 5 5 5 ]\n[ 5 5 5 ]\n");
    EXPECT(A+B,"(1:3),(1:3) \n[ 18 18 18 ]\n[ 18 18 18 ]\n[ 18 18 18 ]\n");
    EXPECT(A-B,"(1:3),(1:3) \n[ 8 8 8 ]\n[ 8 8 8 ]\n[ 8 8 8 ]\n");
    EXPECT(Dot(A,B),"585");
    EXPECT(A+2.0,"(1:3),(1:3) \n[ 15 15 15 ]\n[ 15 15 15 ]\n[ 15 15 15 ]\n");
    EXPECT(A-2.0,"(1:3),(1:3) \n[ 11 11 11 ]\n[ 11 11 11 ]\n[ 11 11 11 ]\n");
    EXPECT(A*2.0,"(1:3),(1:3) \n[ 26 26 26 ]\n[ 26 26 26 ]\n[ 26 26 26 ]\n");
    EXPECT(A/2.0,"(1:3),(1:3) \n[ 6.5 6.5 6.5 ]\n[ 6.5 6.5 6.5 ]\n[ 6.5 6.5 6.5 ]\n");
    EXPECT(A+=B,"(1:3),(1:3) \n[ 18 18 18 ]\n[ 18 18 18 ]\n[ 18 18 18 ]\n");
    EXPECT(A-=B,"(1:3),(1:3) \n[ 13 13 13 ]\n[ 13 13 13 ]\n[ 13 13 13 ]\n");
    EXPECT(A+=2.0,"(1:3),(1:3) \n[ 15 15 15 ]\n[ 15 15 15 ]\n[ 15 15 15 ]\n");
    EXPECT(A-=2.0,"(1:3),(1:3) \n[ 13 13 13 ]\n[ 13 13 13 ]\n[ 13 13 13 ]\n");
    EXPECT(A*=2.0,"(1:3),(1:3) \n[ 26 26 26 ]\n[ 26 26 26 ]\n[ 26 26 26 ]\n");
    EXPECT(A/=2.0,"(1:3),(1:3) \n[ 13 13 13 ]\n[ 13 13 13 ]\n[ 13 13 13 ]\n");
    EXPECT(DirectMultiply(A,B),"(1:3),(1:3) \n[ 65 65 65 ]\n[ 65 65 65 ]\n[ 65 65 65 ]\n");
    EXPECT(DirectDivide(A,B),"(1:3),(1:3) \n[ 2.6 2.6 2.6 ]\n[ 2.6 2.6 2.6 ]\n[ 2.6 2.6 2.6 ]\n");
    EXPECT(A==A,"1");
    EXPECT(A==B,"0");
    EXPECT(A!=A,"0");
    EXPECT(A!=B,"1");
    EXPECT(Sum(A),"117");
    EXPECT(Sum(pow2(A)),"1521");
    //
    //  Verify that mixing types also works.
    //
    EXPECT(A+2,"(1:3),(1:3) \n[ 15 15 15 ]\n[ 15 15 15 ]\n[ 15 15 15 ]\n");
    EXPECT(A-2,"(1:3),(1:3) \n[ 11 11 11 ]\n[ 11 11 11 ]\n[ 11 11 11 ]\n");
    EXPECT(A*2,"(1:3),(1:3) \n[ 26 26 26 ]\n[ 26 26 26 ]\n[ 26 26 26 ]\n");
    EXPECT(A/2,"(1:3),(1:3) \n[ 6.5 6.5 6.5 ]\n[ 6.5 6.5 6.5 ]\n[ 6.5 6.5 6.5 ]\n");
    EXPECT(A+=2,"(1:3),(1:3) \n[ 15 15 15 ]\n[ 15 15 15 ]\n[ 15 15 15 ]\n");
    EXPECT(A-=2,"(1:3),(1:3) \n[ 13 13 13 ]\n[ 13 13 13 ]\n[ 13 13 13 ]\n");
    EXPECT(A*=2,"(1:3),(1:3) \n[ 26 26 26 ]\n[ 26 26 26 ]\n[ 26 26 26 ]\n");
    EXPECT(A/=2,"(1:3),(1:3) \n[ 13 13 13 ]\n[ 13 13 13 ]\n[ 13 13 13 ]\n");

    EXPECT3(A,FillRandom(A,100.0),"(1:3),(1:3) \n[ * * * ]\n[ * * * ]\n[ * * * ]\n");
    EXPECT(Max(A)<=100,"1");
    EXPECT(Max(A)>=0,"1");
    return pass;
}

bool TestM4()
{
    bool pass=true;

    MatrixType Ran(100,10);
    FillRandom(Ran,100.0);
    double min=Ran(1,1), max=Ran(1,1);
    for (int i=1; i<=100; i++)
        for (int j=1; j<=10; j++)
        {
            min=Min(min,Ran(i,j));
            max=Max(max,Ran(i,j));
        }
    EXPECT(Min(Ran),min);
    EXPECT(Max(Ran),max);
    return pass;
}


//
//  Matrix algebra.
//
bool TestM5()
{
    bool pass=true;

    MatrixType A(3,3),C(3,3);
    FillLinear(A,0.0,8.0);
    MatrixType B(Transpose(A));
    EXPECT(A,"(1:3),(1:3) \n[ 0 3 6 ]\n[ 1 4 7 ]\n[ 2 5 8 ]\n");
    EXPECT(B,"(1:3),(1:3) \n[ 0 1 2 ]\n[ 3 4 5 ]\n[ 6 7 8 ]\n");
    EXPECT3(C,C=A*B,"(1:3),(1:3) \n[ 45 54 63 ]\n[ 54 66 78 ]\n[ 63 78 93 ]\n");
    EXPECT3(C,C=-A*B,"(1:3),(1:3) \n[ -45 -54 -63 ]\n[ -54 -66 -78 ]\n[ -63 -78 -93 ]\n");
    EXPECT3(C,C=A*-B,"(1:3),(1:3) \n[ -45 -54 -63 ]\n[ -54 -66 -78 ]\n[ -63 -78 -93 ]\n");
    EXPECT3(C,C=-A*-B,"(1:3),(1:3) \n[ 45 54 63 ]\n[ 54 66 78 ]\n[ 63 78 93 ]\n");
    EXPECT3(C,C=A*B-B*A,"(1:3),(1:3) \n[ 40 40 40 ]\n[ 40 16 -8 ]\n[ 40 -8 -56 ]\n");
    VectorType V(3),VC(3);
    FillLinear(V,0.0,2.0);
    EXPECT(V,"(1:3){ 0 1 2 }");
    EXPECT3(VC,VC=A*V,"(1:3){ 15 18 21 }");
    EXPECT3(VC,VC=V*A,"(1:3){ 5 14 23 }");
    EXPECT(V*A*V,60);
    EXPECT3(VC,VC=-A*V,"(1:3){ -15 -18 -21 }");
    EXPECT3(VC,VC=V*-A,"(1:3){ -5 -14 -23 }");
    EXPECT3(VC,VC=A*-V,"(1:3){ -15 -18 -21 }");
    EXPECT3(VC,VC=-V*A,"(1:3){ -5 -14 -23 }");
    EXPECT3(VC,VC=-A*-V,"(1:3){ 15 18 21 }");
    EXPECT3(VC,VC=-V*-A,"(1:3){ 5 14 23 }");

    EXPECT(A,"(1:3),(1:3) \n[ 0 3 6 ]\n[ 1 4 7 ]\n[ 2 5 8 ]\n");
    EXPECT(B,"(1:3),(1:3) \n[ 0 1 2 ]\n[ 3 4 5 ]\n[ 6 7 8 ]\n");
    EXPECT(A+=~B,"(1:3),(1:3) \n[ 0 6 12 ]\n[ 2 8 14 ]\n[ 4 10 16 ]\n");
    EXPECT(A-=~B,"(1:3),(1:3) \n[ 0 3 6 ]\n[ 1 4 7 ]\n[ 2 5 8 ]\n");
    EXPECT(A+=OuterProduct(V,V),"(1:3),(1:3) \n[ 0 3 6 ]\n[ 1 5 9 ]\n[ 2 7 12 ]\n");
    EXPECT(A-=OuterProduct(V,V),"(1:3),(1:3) \n[ 0 3 6 ]\n[ 1 4 7 ]\n[ 2 5 8 ]\n");
    EXPECT(A+=OuterProduct(V+V,V),"(1:3),(1:3) \n[ 0 3 6 ]\n[ 1 6 11 ]\n[ 2 9 16 ]\n");
    EXPECT(A-=OuterProduct(V+V,V),"(1:3),(1:3) \n[ 0 3 6 ]\n[ 1 4 7 ]\n[ 2 5 8 ]\n");
    EXPECT(A+=OuterProduct(V,V+V),"(1:3),(1:3) \n[ 0 3 6 ]\n[ 1 6 11 ]\n[ 2 9 16 ]\n");
    EXPECT(A-=OuterProduct(V,V+V),"(1:3),(1:3) \n[ 0 3 6 ]\n[ 1 4 7 ]\n[ 2 5 8 ]\n");
    EXPECT(A+=OuterProduct(2.0*V,V+V),"(1:3),(1:3) \n[ 0 3 6 ]\n[ 1 8 15 ]\n[ 2 13 24 ]\n");
    EXPECT(A-=OuterProduct(2.0*V,V+V),"(1:3),(1:3) \n[ 0 3 6 ]\n[ 1 4 7 ]\n[ 2 5 8 ]\n");
    EXPECT(A+=2.0*OuterProduct(V,V),"(1:3),(1:3) \n[ 0 3 6 ]\n[ 1 6 11 ]\n[ 2 9 16 ]\n");
    EXPECT(A-=OuterProduct(V,V)*2.0,"(1:3),(1:3) \n[ 0 3 6 ]\n[ 1 4 7 ]\n[ 2 5 8 ]\n");

    EXPECT3(A,Unit(A),"(1:3),(1:3) \n[ 1 0 0 ]\n[ 0 1 0 ]\n[ 0 0 1 ]\n");
    EXPECT(IsSymmetric(A),true);
    EXPECT(IsSymmetric(B),false);
    EXPECT(MakeSymmetric(B),4);
    EXPECT(B,"(1:3),(1:3) \n[ 0 2 4 ]\n[ 2 4 6 ]\n[ 4 6 8 ]\n");

    return pass;
}


//-----------------------------------------------------
//
//  IO tests
//
bool TestM6()
{
    bool pass=true;

    StreamableObject::SetToAscii();
    {
        MatrixType A(-10,20,-5,50),B;
        FillRandom(A,100.0);
        {
            std::ofstream file("temp.dat");
            file << A;
        }
        StreamableObject::SetToPretty();
        {
            std::ifstream file("temp.dat");
            file >> B;
        }
        EXPECT1(B,A,"file >> B in ascii mode");
    }

    StreamableObject::SetToBinary();
    {
        MatrixType A(-10,20,-5,5),B;
        FillRandom(A,100.0);
        {
            std::ofstream file("tempb.dat",std::ios::out|std::ios::binary);
            file << A;
        }
        StreamableObject::SetToPretty();
        {
            std::ifstream file("tempb.dat",std::ios::in|std::ios::binary);
            file >> B;
        }
        EXPECT1(B,A,"file >> B in binary mode");
    }
    return pass;
}

bool TestM7()
{
    bool pass=true;

    int N=200;
    MatrixType A(N,N),B(N,N),C(N,N);
    VectorType V(N);
    FillRandom(A,100.0);
    FillRandom(B,100.0);
    FillRandom(V,100.0);

    StopWatch s;
    s.Start();
    for (int i=1; i<=10; i++) C=A*B;
    s.Stop();
    double Flops=10*N*N*N*2/s.GetTime();
    std::cout << "Matrix multiply C=A*B Exp template " << Flops*1e-6 << " Mips" << std::endl;
    s.Start();
    for (int i=1; i<=10; i++) mmul(C,A,B);
    s.Stop();
    Flops=10*N*N*N*2/s.GetTime();
    std::cout << "Matrix multiply C=A*B Hand coaded  " << Flops*1e-6 << " Mips" << std::endl;

    N=1000;
    int N1=10;
    A.SetLimits(N,N);
    V.SetLimits(N);
    FillRandom(A,10.0);
    FillRandom(V,10.0);
    double d1, d2=0;//Warning bug
    s.Start();
    for (int i=1; i<=N1; i++) d1=vmul1(V,A,V);
    s.Stop();

    double MFlops=(N1*1e-6)*2*(N*N+N)/s.GetTime();
    std::cout << "Matrix multiply V*M*V Exp template " << MFlops << " Mips" << std::endl;
    s.Start();
    for (int i=1; i<=N1; i++) d2=vmul(V,A,V);
    s.Stop();
    assert(fabs(d1-d2)<1e-6);
    MFlops=(N1*1e-6)*2*(N*N+N)/s.GetTime();
    std::cout << "Matrix multiply V*M*V Hand coaded  " << MFlops << " Mips" << std::endl;
    return pass;
}

bool TestM8()
{

    int N=200;
    double rmax=10.0;
    MatrixType A(N,N),B(N,N),C(N,N);
    VectorType V(N);
    FillRandom(A,rmax);
    FillRandom(B,rmax);
    FillRandom(V,rmax);

    StopWatch s;
    s.Start();
    for (int i=1; i<=10; i++) C=A*B;
    s.Stop();
    double Flops=10*N*N*N*2/s.GetTime();
    std::cout << "Matrix multiply C=A*B Exp template " << Flops*1e-6 << " MFlops" << std::endl;
    s.Start();
    for (int i=1; i<=10; i++) mmul(C,A,B);
    s.Stop();
    Flops=10*N*N*N*2/s.GetTime();
    std::cout << "Matrix multiply C=A*B Hand coaded  " << Flops*1e-6 << " MFlops" << std::endl;

    N=1000;
    int N1=50;
    A.SetLimits(N,N);
    V.SetLimits(N);
    FillRandom(A,rmax);
    FillRandom(V,rmax);
    double d;
    s.Start();
    for (int i=1; i<=N1; i++) d=V*A*V;
    s.Stop();
    std::cout << d << std::endl;
    double MFlops=(N1*1e-6)*2*(N*N+N)/s.GetTime();
    std::cout << "Matrix std::vector multiply V*M*V Exp template " << MFlops << " MFlops" << std::endl;
    s.Start();
    for (int i=1; i<=N1; i++) d=vmul(V,A,V);
    s.Stop();
    std::cout << d << std::endl;
    MFlops=(N1*1e-6)*2*(N*N+N)/s.GetTime();
    std::cout << "Matrix std::vector multiply V*M*V Hand coaded  " << MFlops << " MFlops" << std::endl;


    return true;

}

bool TestM9()
{
    bool pass=true;
    {
        SMatrix<double> SRan(10,10),SOpRan(10,10);
        MatrixType Ran(10,10),OpRan(10,10),U(10,10),A(10,10);
        FillPos(Ran);
        FillPos(SRan);
        Unit(U);

        OpRan=Ran*InvertSymmetric(Ran);
        EXPECT(Sum(abs(OpRan-U))<1e-10,true);

        OpRan=SRan*InvertSymmetric(SRan);
        EXPECT(Sum(abs(OpRan-U))<1e-10,true);

        FillRandomPositive(Ran,10.0);
        Ran=(Ran+Matrix<double>(~Ran))/2.0;
        OpRan=Ran*InvertSquare(Ran);
        EXPECT(Sum(abs(OpRan-U))<1e-5,true);

        Vector<double> X(10),V(10);
        Matrix<double> Temp;
        Matrix<double> B(10,10),BTemp(10,10);
        FillRandom(A,10.0);
        FillRandom(V,10.0);
        Temp=A;
        X=V;
        LUSolver<double> lu(A);
        lu.SolveFor(V);
        EXPECT(Sum(abs(Temp*V-X))<1e-10,true);

//    FillRandom(B,10.0);
        Unit(B);
        BTemp=B;
        lu.SolveFor(B);
        EXPECT(Sum(abs(Matrix<double>(Temp*B)-BTemp))<1e-10,true);

        for (int i=1; i<=10; i++)
        {
            V=B.GetColumn(i);
            X=V;
            lu.SolveFor(V);
            EXPECT(Sum(abs(Temp*V-X))<1e-10,true);
        }

        FillRandom(A,10.0);
        FillRandom(V,10.0);
        Temp=A;
        X=V;
        SVSolver<double> sv(A);
        V=sv.SolveFor(X);
        EXPECT(Sum(abs(Temp*V-X))<1e-10,true);

        FillRandom(A,10.0);
        FillRandom(B,10.0);
        Temp=A;
        BTemp=B;
        GaussJordanSolver(A,B);
        EXPECT(Sum(abs(Matrix<double>(Temp*B)-BTemp))<1e-10,true);

        FillPos(Ran);
        OpRan=Ran;
        Cholsky(OpRan);
        EXPECT(Sum(abs(Matrix<double>(OpRan*~OpRan)-Ran))<1e-10,true);

        OpRan=Ran;
        Vector<double> EigenValues=Diagonalize(OpRan);
        Matrix<double> Diag=~OpRan*Ran*OpRan;
        EXPECT(Sum(abs(Diag.GetDiagonal()-EigenValues))<1e-10,true);
    }
    {
        int N=300;
        Matrix<double> Ran(N,N);
        FillRandom(Ran,10.0);
        MakeSymmetric(Ran);

        StopWatch s;
        s.Start();
        Vector<double>  EigenValues=Diagonalize(Ran);
        s.Stop();
        std::cout << s.GetTime() << " seconds to diagonalize a " << N << "X" << N << " real matrix." << std::endl;
    }
//-----------------------------------------------------
//
//  IO tests
//
    StreamableObject::SetToAscii();
    {
        MatrixType A(-10,20,-5,50),B;
        FillRandom(A,100.);
        {
            std::ofstream file("temp.dat");
            file << A;
        }
        StreamableObject::SetToPretty();
        {
            std::ifstream file("temp.dat");
            file >> B;
        }
        EXPECT1(B,A,"file >> B in ascii mode");
    }

    StreamableObject::SetToBinary();
    {
        MatrixType A(-10,20,-5,5),B;
        FillRandom(A,100.);
        {
            std::ofstream file("tempb.dat",std::ios::out|std::ios::binary);
            file << A;
        }
        StreamableObject::SetToPretty();
        {
            std::ifstream file("tempb.dat",std::ios::in|std::ios::binary);
            file >> B;
        }
        EXPECT1(B,A,"file >> B in binary mode");
    }
    return pass;
}
