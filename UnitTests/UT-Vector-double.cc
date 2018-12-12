// File: UT-Vector-double.cc  Unit test the Vector class for double data types.

// Copyright (1994-2005), Jan N. Reimers

#include "oml/UnitTests/UnitTest.H"
#include "oml/vector.h"
#include "oml/vector_io.h"
#include "oml/shellsort.h"
#include "oml/heapsort.h"
#include "oml/indexsort.h"
#include "oml/random.h"
#include <cmath>
#include <fstream>

  bool IsASorted(const Vector<double>& a)
  {
    bool ret=true;
    for (index_t i=a.GetLow()+1; i<=a.GetHigh(); i++) ret=ret&&(a(i-1)<=a(i));
    return ret;
  }

  bool IsDSorted(const Vector<double>& a)
  {
    bool ret=true;
    for (index_t i=a.GetLow()+1; i<=a.GetHigh(); i++) ret=ret&&(a(i-1)>=a(i));
    return ret;
  }

template <class T> T MaxAbs(const Vector<T>& v)
{
  return Max(abs(v));
}

int TestVectorDouble()
{
  const char* Class="Vector<double>";
  bool pass=true;
//
//  Start message.
//

  StartClass(Class);
  StreamableObject::SetOutputMode(StreamableObject::pretty);
	std::cout << "Vector<double> foot print=" << sizeof(Vector<double>) << " bytes" << std::endl;
  {
    Vector<double> A(-1,8);
    EXPECT3(A,FillPower(A,1.0,10.0),"(-1:8){ 1 1.29155 1.6681 2.15443 2.78256 3.59381 4.64159 5.99484 7.74264 10 }");
    EXPECT3(A,FillRandom(A,100.0),"(-1:8){ * * * * * * * * * * }");
    EXPECT(Max(A)<= 100,true);
    EXPECT(Min(A)>=-100,true);
  }
  {
    Vector<double> Ran(-1,1000);
    FillRandom(Ran,100.0);
    double min=Ran(-1), max=Ran(-1);
    for (int i=-1;i<1000;i++)
    {
      min=Min(min,Ran(i));
      max=Max(max,Ran(i));
    }
    EXPECT(Min(Ran),min);
    EXPECT(Max(Ran),max);
  }
  {  //Mixing double and ints.
    Vector<double> A(6);
    Fill(A,13.0);
    EXPECT(A+2,"(1:6){ 15 15 15 15 15 15 }");
    EXPECT(A-2,"(1:6){ 11 11 11 11 11 11 }");
    EXPECT(A*2,"(1:6){ 26 26 26 26 26 26 }");
    EXPECT(A/2,"(1:6){ 6.5 6.5 6.5 6.5 6.5 6.5 }");
    EXPECT(2+A,"(1:6){ 15 15 15 15 15 15 }");
    EXPECT(2-A,"(1:6){ -11 -11 -11 -11 -11 -11 }");
    EXPECT(2*A,"(1:6){ 26 26 26 26 26 26 }");
    EXPECT(26/A,"(1:6){ 2 2 2 2 2 2 }");
  }
  {
    Vector<double> A(10),B(10);
    double Pi=M_PI;
    Fill(A,Pi);
    Fill(B,0.25);
    EXPECT(MaxAbs(Vector<double>(1.0/A   -(1.0/ Pi)))<eps,true);
    EXPECT(MaxAbs(Vector<double>(sin(A)  - sin (Pi)))<eps,true);
    EXPECT(MaxAbs(Vector<double>(cos(A)  - cos (Pi)))<eps,true);
    EXPECT(MaxAbs(Vector<double>(tan(A)  - tan (Pi)))<eps,true);
    EXPECT(MaxAbs(Vector<double>(asin(B) -asin (0.25)))<eps,true);
    EXPECT(MaxAbs(Vector<double>(acos(B) -acos (0.25)))<eps,true);
    EXPECT(MaxAbs(Vector<double>(atan(A) -atan (Pi)))<eps,true);
    EXPECT(MaxAbs(Vector<double>(sinh(A) -sinh (Pi)))<eps,true);
    EXPECT(MaxAbs(Vector<double>(cosh(A) -cosh (Pi)))<eps,true);
    EXPECT(MaxAbs(Vector<double>(tanh(A) -tanh (Pi)))<eps,true);
    EXPECT(MaxAbs(Vector<double>(log(A)  -log  (Pi)))<eps,true);
    EXPECT(MaxAbs(Vector<double>(log10(A)-log10(Pi)))<eps,true);
    EXPECT(MaxAbs(Vector<double>(exp(A)  -exp  (Pi)))<eps,true);
    EXPECT(MaxAbs(Vector<double>(sqrt(A) -sqrt (Pi)))<eps,true);
    EXPECT(MaxAbs(Vector<double>(pow2(A) -    Pi*Pi))<eps,true);
    EXPECT(MaxAbs(Vector<double>(pow3(A) - Pi*Pi*Pi))<eps,true);
    EXPECT(MaxAbs(Vector<double>(pow4(A)  -Pi*Pi*Pi*Pi))<100*eps,true);
    //EXPECT(MaxAbs(Vector<double>(pow10(A)  -pow10(Pi)))<1000*eps,true);
    EXPECT(MaxAbs(Vector<double>(pow(A,0.6)-pow   (Pi,0.6)))<eps,true);
    EXPECT(MaxAbs(Vector<double>(pow(A,B)  -pow  (Pi,0.25)))<eps,true);
    EXPECT(MaxAbs(Vector<double>(atan2(A,B)-atan2(Pi,0.25)))<eps,true);
  }
//##########################################################################
//
//  Test the sorting routines
//
  {
    Vector<double> A(1000),B(1000);
    FillRandom(A,100.0);
    B=A;
    EXPECT3(IsASorted(B),AscendingShellSort(B),true);
    B=A;
    EXPECT3(IsDSorted(B),DescendingShellSort(B),true);
    B=A;
    EXPECT3(IsASorted(B),AscendingHeapSort(B),true);
    B=A;
    EXPECT3(IsDSorted(B),DescendingHeapSort(B),true);
    B=A;
    EXPECT3(IsASorted(B),B.ReIndex(MakeAscendingIndex(A)),true);
    B=A;
    EXPECT3(IsDSorted(B),B.ReIndex(MakeDescendingIndex(A)),true);

  }
  {
    Vector<double> A(10);
    Fill(A,M_PI);
    EXPECT(fabs(!A-sqrt(10.0*M_PI*M_PI))<1e-15,true);
    EXPECT3(A*A,Normalize(A),1.0);
  }

//-----------------------------------------------------
//
//  IO tests
//
  StreamableObject::SetToAscii();
  {
    Vector<double> A(-10,1000),B;
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
    Vector<double> A(-10,1000),B;
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
  return pass ? 0 : -1;
}



