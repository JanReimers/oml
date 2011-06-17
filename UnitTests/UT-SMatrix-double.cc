// File: UT-SMatrix-double.cc  Unit test the SMatrix class for double data types.

// Copyright (1994-2005), Jan N. Reimers

#include "oml/UnitTests/UnitTest.H"
#include "oml/vector_io.h"
#include "oml/smatrix.h"
#include "oml/matrix.h"
#include "oml/minmax.h"
#include "oml/random.h"
#include "oml/stopw.h"
#include <fstream>

typedef SMatrix<double> MatrixType;
typedef Vector <double> VectorType;


void mmul(Matrix<double> & C,const MatrixType& A, const MatrixType&B);
double vmul(const VectorType V1,const MatrixType& A, const VectorType& V2);
double vmul1(const VectorType V1,const MatrixType& A, const VectorType& V2);


int TestSMatrixDouble()
{
  const char* Class="SymMatrix<double>";
  bool pass=true;
//
//  Start message.
//

  StartClass(Class);
  StreamableObject::SetToPretty();

//###########################################################################
//
//  Constructors
//
  {
    EXPECT(MatrixType(),"(1:0),(1:0) ");
    EXPECT(MatrixType(3,3),"(1:3),(1:3) \n[ * * * ]\n[ * * ]\n[ * ]\n");
    EXPECT(MatrixType(-1,3,-1,3),
	   "(-1:3),(-1:3) \n[ * * * * * ]\n[ * * * * ]\n[ * * * ]\n[ * * ]\n[ * ]\n");
    EXPECT(MatrixType(VecLimits(-1,3),VecLimits(-1,3)),
	   "(-1:3),(-1:3) \n[ * * * * * ]\n[ * * * * ]\n[ * * * ]\n[ * * ]\n[ * ]\n");
    EXPECT(MatrixType(MatLimits(-1,3,-1,3)),
	   "(-1:3),(-1:3) \n[ * * * * * ]\n[ * * * * ]\n[ * * * ]\n[ * * ]\n[ * ]\n");
    MatrixType A1(3,3),A2(3,3);
    Fill(A1,3.0);
    EXPECT3(temp,MatrixType temp(A1),"(1:3),(1:3) \n[ 3 3 3 ]\n[ 3 3 ]\n[ 3 ]\n");
    EXPECT(A2=A1,"(1:3),(1:3) \n[ 3 3 3 ]\n[ 3 3 ]\n[ 3 ]\n");

    EXPECT2(A1,SetLimits(0,4,0,4,true),
		"(0:4),(0:4) \n[ * * * * * ]\n[ 3 3 3 * ]\n[ 3 3 * ]\n[ 3 * ]\n[ * ]\n");
    EXPECT2(A1,SetLimits(0,4,0,4),
		"(0:4),(0:4) \n[ * * * * * ]\n[ * * * * ]\n[ * * * ]\n[ * * ]\n[ * ]\n");

    EXPECT3(A2,FillLinear(A2,1.0,6.0),"(1:3),(1:3) \n[ 1 2 3 ]\n[ 4 5 ]\n[ 6 ]\n");
    EXPECT(A2(2,3),5);
		A2(3,2)=9;
		EXPECT3(A2,A2(3,2)=9,"(1:3),(1:3) \n[ 1 2 3 ]\n[ 4 9 ]\n[ 6 ]\n");
    EXPECT(A2.SubMatrix(MatLimits(2,3,2,3)),"(2:3),(2:3) \n[ 4 9 ]\n[ 6 ]\n");
	}
  {
    MatrixType A(2,8,2,8);
    FillLinear(A,1.0,28.0);
//		std::cout << A;
    EXPECT(A.GetDiagonal(),"(2:8){ 1 8 14 19 23 26 28 }");
    EXPECT(A.GetColumn (3),"(2:8){ 2 8 9 10 11 12 13 }");
    EXPECT(A.GetRow    (6),"(2:8){ 5 11 16 20 23 24 25 }");
    EXPECT(Sum(A.GetColumn (3)),65);
    EXPECT(Sum(A.GetRow    (6)),124);
		EXPECT(Sum(A.GetDiagonal()),119);
		EXPECT(A.GetRow(6)*A.GetColumn (3),1308);

  }
  {
    MatrixType A(3,3),B(3,3);
    Fill(A,13.0);
    Fill(B, 5.0);

    EXPECT(A,"(1:3),(1:3) \n[ 13 13 13 ]\n[ 13 13 ]\n[ 13 ]\n");
    EXPECT(B,"(1:3),(1:3) \n[ 5 5 5 ]\n[ 5 5 ]\n[ 5 ]\n");
    EXPECT(A+B,"(1:3),(1:3) \n[ 18 18 18 ]\n[ 18 18 ]\n[ 18 ]\n");
    EXPECT(A-B,"(1:3),(1:3) \n[ 8 8 8 ]\n[ 8 8 ]\n[ 8 ]\n");
    EXPECT(Dot(A,B),"585");
    EXPECT(A+2.0,"(1:3),(1:3) \n[ 15 15 15 ]\n[ 15 15 ]\n[ 15 ]\n");
    EXPECT(A-2.0,"(1:3),(1:3) \n[ 11 11 11 ]\n[ 11 11 ]\n[ 11 ]\n");
    EXPECT(A*2.0,"(1:3),(1:3) \n[ 26 26 26 ]\n[ 26 26 ]\n[ 26 ]\n");
    EXPECT(A/2.0,"(1:3),(1:3) \n[ 6.5 6.5 6.5 ]\n[ 6.5 6.5 ]\n[ 6.5 ]\n");
    EXPECT(A+=B,"(1:3),(1:3) \n[ 18 18 18 ]\n[ 18 18 ]\n[ 18 ]\n");
    EXPECT(A-=B,"(1:3),(1:3) \n[ 13 13 13 ]\n[ 13 13 ]\n[ 13 ]\n");
    EXPECT(A+=2.0,"(1:3),(1:3) \n[ 15 15 15 ]\n[ 15 15 ]\n[ 15 ]\n");
    EXPECT(A-=2.0,"(1:3),(1:3) \n[ 13 13 13 ]\n[ 13 13 ]\n[ 13 ]\n");
    EXPECT(A*=2.0,"(1:3),(1:3) \n[ 26 26 26 ]\n[ 26 26 ]\n[ 26 ]\n");
    EXPECT(A/=2.0,"(1:3),(1:3) \n[ 13 13 13 ]\n[ 13 13 ]\n[ 13 ]\n");
    EXPECT(DirectMultiply(A,B),"(1:3),(1:3) \n[ 65 65 65 ]\n[ 65 65 ]\n[ 65 ]\n");
    EXPECT(DirectDivide(A,B),"(1:3),(1:3) \n[ 2.6 2.6 2.6 ]\n[ 2.6 2.6 ]\n[ 2.6 ]\n");
    EXPECT(A==A,"1");
    EXPECT(A==B,"0");
    EXPECT(A!=A,"0");
    EXPECT(A!=B,"1");
    EXPECT(Sum(A),"117");
    EXPECT(Sum(pow2(A)),"1521");

    EXPECT3(A,FillRandom(A,100.0),"(1:3),(1:3) \n[ * * * ]\n[ * * ]\n[ * ]\n");
    EXPECT(Max(A)<=100,"1");
    EXPECT(Max(A)>=0,"1");

    {
      MatrixType Ran(100,100);
      FillRandom(Ran,100.0);
      double min=Ran(1,1), max=Ran(1,1);
      for (int i=1;i<=100;i++)
      for (int j=1;j<=100;j++)
      {
	min=Min(min,Ran(i,j));
	max=Max(max,Ran(i,j));
      }
      EXPECT(Min(Ran),min);
      EXPECT(Max(Ran),max);
    }
  }
//
//  Matrix algebra.
//
  {
    MatrixType A(3,3),B(3,3);
    Matrix<double> C(3,3);
    FillLinear(A,0.0,5.0);
    FillLinear(B,5.0,0.0);

    EXPECT(A,"(1:3),(1:3) \n[ 0 1 2 ]\n[ 3 4 ]\n[ 5 ]\n");
    EXPECT(B,"(1:3),(1:3) \n[ 5 4 3 ]\n[ 2 1 ]\n[ 0 ]\n");
    EXPECT3(C,C=A*B,"(1:3),(1:3) \n[ 10 4 1 ]\n[ 29 14 6 ]\n[ 41 21 10 ]\n");
    EXPECT3(C,C=-A*B,"(1:3),(1:3) \n[ -10 -4 -1 ]\n[ -29 -14 -6 ]\n[ -41 -21 -10 ]\n");
    EXPECT3(C,C=A*-B,"(1:3),(1:3) \n[ -10 -4 -1 ]\n[ -29 -14 -6 ]\n[ -41 -21 -10 ]\n");
    EXPECT3(C,C=-A*-B,"(1:3),(1:3) \n[ 10 4 1 ]\n[ 29 14 6 ]\n[ 41 21 10 ]\n");
    EXPECT3(C,C=A*B-B*A,"(1:3),(1:3) \n[ 0 -25 -40 ]\n[ 25 0 -15 ]\n[ 40 15 0 ]\n");
    VectorType V(3),VC(3);
    FillLinear(V,0.0,2.0);
    EXPECT(V,"(1:3){ 0 1 2 }");
    EXPECT3(VC,VC=A*V,"(1:3){ 5 11 14 }");
    EXPECT3(VC,VC=V*A,"(1:3){ 5 11 14 }");
    EXPECT(V*A*V,39);
    EXPECT3(VC,VC=-A*V,"(1:3){ -5 -11 -14 }");
    EXPECT3(VC,VC=V*-A,"(1:3){ -5 -11 -14 }");
    EXPECT3(VC,VC=A*-V,"(1:3){ -5 -11 -14 }");
    EXPECT3(VC,VC=-V*A,"(1:3){ -5 -11 -14 }");
    EXPECT3(VC,VC=-A*-V,"(1:3){ 5 11 14 }");
    EXPECT3(VC,VC=-V*-A,"(1:3){ 5 11 14 }");

		EXPECT(A,"(1:3),(1:3) \n[ 0 1 2 ]\n[ 3 4 ]\n[ 5 ]\n");
    EXPECT(B,"(1:3),(1:3) \n[ 5 4 3 ]\n[ 2 1 ]\n[ 0 ]\n");
    EXPECT(A+=OuterProduct(V),"(1:3),(1:3) \n[ 0 1 2 ]\n[ 4 6 ]\n[ 9 ]\n");
    EXPECT(A-=OuterProduct(V),"(1:3),(1:3) \n[ 0 1 2 ]\n[ 3 4 ]\n[ 5 ]\n");
    EXPECT(A+=OuterProduct(V+V),"(1:3),(1:3) \n[ 0 1 2 ]\n[ 7 12 ]\n[ 21 ]\n");
    EXPECT(A-=OuterProduct(V+V),"(1:3),(1:3) \n[ 0 1 2 ]\n[ 3 4 ]\n[ 5 ]\n");
    EXPECT(A+=2.0*OuterProduct(V),"(1:3),(1:3) \n[ 0 1 2 ]\n[ 5 8 ]\n[ 13 ]\n");
    EXPECT(A-=OuterProduct(V)*2.0,"(1:3),(1:3) \n[ 0 1 2 ]\n[ 3 4 ]\n[ 5 ]\n");

    EXPECT3(A,Unit(A),"(1:3),(1:3) \n[ 1 0 0 ]\n[ 1 0 ]\n[ 1 ]\n");

  }
//-----------------------------------------------------
//
//  IO tests
//
  StreamableObject::SetToAscii();
  {
    MatrixType A(-10,20,-10,20),B;
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
    MatrixType A(-10,20,-10,20),B;
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




  {
    int N=200;
    MatrixType A(N,N),B(N,N);
    Matrix<double> C(N,N);
    VectorType V(N);
    FillRandom(A,100.0);
    FillRandom(B,100.0);
    FillRandom(V,100.0);

    StopWatch s;
    s.Start();
    for (int i=1;i<=10;i++) C=A*B;
    s.Stop();
    double Flops=10*N*N*N*2/s.GetTime();
    std::cout << "Symmetric Matrix multiply C=A*B Exp template " << Flops*1e-6 << " Mips" << std::endl;
    s.Start();
    for (int i=1;i<=10;i++) mmul(C,A,B);
    s.Stop();
    Flops=10*N*N*N*2/s.GetTime();
    std::cout << "Symmetric Matrix multiply C=A*B Hand coaded  " << Flops*1e-6 << " Mips" << std::endl;

    N=1000;
    int N1=10;
    A.SetLimits(N,N);
    V.SetLimits(N);
    FillRandom(A,10.0);
    FillRandom(V,10.0);
    double d1,d2=0;
    s.Start();
    for (int i=1;i<=N1;i++) d1=vmul1(V,A,V);
    s.Stop();
    double MFlops=(N1*1e-6)*2*(N*N+N)/s.GetTime();
    std::cout << "Symmetric Matrix multiply V*M*V Exp template " << MFlops << " Mips" << std::endl;
    s.Start();
    for (int i=1;i<=N1;i++) d2=vmul(V,A,V);
    s.Stop();
    assert(fabs(d1-d2)<1e-6);
    MFlops=(N1*1e-6)*2*(N*N+N)/s.GetTime();
    std::cout << "Symmetric Matrix multiply V*M*V Hand coaded  " << MFlops << " Mips" << std::endl;
  }
  return pass ? 0 : -1;
}



void mmul(Matrix<double>& C,const MatrixType& A, const MatrixType&B)
{
  int N=A.GetLimits().Row.High;
  Matrix<double>::Subscriptor s(C);
  for (int i=1;i<=N;i++)
    for (int j=1;j<=N;j++)
    {
      double t=0;
      for (int k=1;k<=N;k++) t+=A(i,k)*B(k,j);
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
  for (int i=1;i<=N;i++)
  {
    double t=0;
    for (int j=1;j<=N;j++) t+=A(i,j)*V2(j);
    ret+=V1(i)*t;
  }
  return ret;
}


