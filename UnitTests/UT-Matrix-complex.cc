// File: UT-Matrix-complex.cc  Unit test the Matrix class for complex data types.

// Copyright (1994-2003), Jan N. Reimers

#include "oml/UnitTests/UnitTest.H"
#include "oml/matrix.h"
#include "oml/smatrix.h"
#include "oml/random.h"
#include "oml/stopw.h"
#include "oml/cnumeric.h"
#include <iomanip>
#include <cmath>
#include <complex>
#include <fstream>

int main()
{
//
//  Start message.
//
  StartClass("Matrix<std::complex<T> >");
  bool pass=true;

  StreamableObject::SetToPretty();
  std::cout.setf(std::ios::fixed,std::ios::floatfield);

  typedef Matrix<std::complex<double> > MatrixType;
  typedef SMatrix<std::complex<double> > HMatrixType;
  
  {
    MatrixType A(3,3);
    HMatrixType B(3,3);
    
    for (int i=1;i<=3;i++)
      for (int j=1;j<=3;j++)
	A(i,j)=std::complex<double>(i+j,i-j);
    
    for (int i=1;i<=3;i++)
      for (int j=i;j<=3;j++)
	B(i,j)=std::complex<double>(i+j,i-j);    

    EXPECT(A,"(1:3),(1:3) \n[ (2,0) (3,-1) (4,-2) ]\n[ (3,1) (4,0) (5,-1) ]\n[ (4,2) (5,1) (6,0) ]\n");
    EXPECT(Sum(A),std::complex<double>(36,0));
    EXPECT(Dot(A,A),std::complex<double>(168,0));    
    EXPECT(Dot(A+A,A),std::complex<double>(336,0));    
    EXPECT(Dot(A,A+A),std::complex<double>(336,0));    
    EXPECT(Dot(A+A,A+A),std::complex<double>(672,0));
    
    EXPECT(B,"(1:3),(1:3) \n[ (2,0) (3,-1) (4,-2) ]\n[ (4,0) (5,-1) ]\n[ (6,0) ]\n");
    const HMatrixType& Br(B);
    EXPECT(Br(1,3)==conj(Br(3,1)),true);
    EXPECT(Sum(B),std::complex<double>(36,0));
    EXPECT(Dot(B,B),std::complex<double>(168,0));    
    EXPECT(Dot(B+B,B),std::complex<double>(336,0));    
    EXPECT(Dot(B,B+B),std::complex<double>(336,0));    
    EXPECT(Dot(B+B,B+B),std::complex<double>(672,0));  
  
  
    EXPECT(A*A,"(1:3),(1:3) \n[ (34,0) (40,-12) (46,-24) ]\n[ (40,12) (52,0) (64,-12) ]\n[ (46,24) (64,12) (82,0) ]\n");
    EXPECT(A*B,"(1:3),(1:3) \n[ (34,0) (40,-12) (46,-24) ]\n[ (40,12) (52,0) (64,-12) ]\n[ (46,24) (64,12) (82,0) ]\n");
    EXPECT(B*A,"(1:3),(1:3) \n[ (34,0) (40,-12) (46,-24) ]\n[ (40,12) (52,0) (64,-12) ]\n[ (46,24) (64,12) (82,0) ]\n");
    EXPECT(B*B,"(1:3),(1:3) \n[ (34,0) (40,-12) (46,-24) ]\n[ (40,12) (52,0) (64,-12) ]\n[ (46,24) (64,12) (82,0) ]\n");
  }
  
  {
		int N=100;
    HMatrixType Ran(N,N);
    MatrixType EigenVectors(N,N);
    FillRandom(Ran,std::complex<double>(10.0,5.0));
    for (int i=1;i<=N;i++) Ran(i,i)=real(Ran(i,i));
    
    StopWatch s;
    s.Start();
    Vector<double>  EigenValues=Diagonalize(EigenVectors,Ran);
    s.Stop();
		std::cout << s.GetTime() << " to diagonalize a " << N << "X" << N << " complex matrix." << std::endl;
    MatrixType      Diag=conj(~EigenVectors)*EigenVectors;
    EXPECT(Sum(abs(imag(Diag)))<1e-10,true);
    EXPECT(Sum(abs(real(Diag.GetDiagonal())-1.0))<1e-10,true);
    MatrixType      Diag1=conj(~EigenVectors)*MatrixType(Ran*EigenVectors);
    EXPECT(Sum(abs(real(Diag1.GetDiagonal())-EigenValues))<1e-10,true);
    EXPECT(Sum(abs(imag(Diag1)))<1e-10,true);
  }
  
//-----------------------------------------------------
//
//  IO tests
//
  StreamableObject::SetToAscii();
  {
    MatrixType A(-10,20,-5,50),B;
    FillRandom(A,std::complex<double>(100));
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
    FillRandom(A,std::complex<double>(100));
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
 //-----------------------------------------------------
//
//  IO tests
//
  StreamableObject::SetToAscii();
  {
    HMatrixType A(-10,20,-10,20),B;
    FillRandom(A,std::complex<double>(100));
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
    HMatrixType A(-10,20,-10,20),B;
    FillRandom(A,std::complex<double>(100));
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
  return pass ? 0 : -1;
}
