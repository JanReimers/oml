// File: UT-List-double.cc  Unit test the List class for double data types.

// Copyright (1994-2005), Jan N. Reimers

#include "oml/UnitTests/UnitTest.H"
#include "oml/list.h"
#include "oml/list_io.h"
#include "oml/array.h"
#include "oml/array_io.h"
#include "oml/shellsort.h"
#include "oml/heapsort.h"
#include "oml/indexsort.h"
#include "oml/random.h"
#include "oml/minmax.h"
#include <fstream>

  bool IsASorted(const List<double>& a)
  {
    bool ret=true;
    for (index_t i=1; i<a.size(); i++) ret=ret&&(a[i-1]<=a[i]);
    return ret;
  }

  bool IsDSorted(const List<double>& a)
  {
    bool ret=true;
    for (index_t i=1; i<a.size(); i++) ret=ret&&(a[i-1]>=a[i]);
    return ret;
  }

int TestListDouble()
{
  const char* Class="List<double>";
  bool pass=true;
//
//  Start message.
//

  StartClass(Class);
  StreamableObject::SetOutputMode(StreamableObject::pretty);
  std::cout << "List<double> foot print=" << sizeof(List<double>) << " bytes" << std::endl;

  {
    List<double> A;
    for (int i=0;i<10;i++) A.push_back(i);
    EXPECT3(A,FillPower(A,1.0,10.0),"{ 1 1.29155 1.6681 2.15443 2.78256 3.59381 4.64159 5.99484 7.74264 10 }");
    EXPECT3(A,FillRandom(A,100.0),"{ * * * * * * * * * * }");
    EXPECT(Max(A)<=100,true);
    EXPECT(Max(A)>=0,true);
  }
  {
    List<double> Ran;
    for (int i=0;i<1000;i++) Ran.push_back(i);
    FillRandom(Ran,100.0);
    double min=Ran[0], max=Ran[0];
    for (int i=0;i<1000;i++)
    {
      min=Min(min,Ran[i]);
      max=Max(max,Ran[i]);
    }
    EXPECT(Min(Ran),min);
    EXPECT(Max(Ran),max);
  }
  {  //Mixing double and ints.
    List<double> A(6);
    Fill(A,13.0);
    EXPECT(A+2,"{ 15 15 15 15 15 15 }");
    EXPECT(A-2,"{ 11 11 11 11 11 11 }");
    EXPECT(A*2,"{ 26 26 26 26 26 26 }");
    EXPECT(A/2,"{ 6.5 6.5 6.5 6.5 6.5 6.5 }");
    EXPECT(2+A,"{ 15 15 15 15 15 15 }");
    EXPECT(2-A,"{ -11 -11 -11 -11 -11 -11 }");
    EXPECT(2*A,"{ 26 26 26 26 26 26 }");
    EXPECT(26/A,"{ 2 2 2 2 2 2 }");
  }
  {
    List<double> A,B;
    double Pi=M_PI;
    for (int i=0;i<10;i++) A.push_back(M_PI);
    for (int i=0;i<10;i++) B.push_back(0.25);
    EXPECT(Max(abs(1.0/A   -(1.0/ Pi)))<eps,true);
    EXPECT(Max(abs(sin(A)  - sin (Pi)))<eps,true);
    EXPECT(Max(abs(cos(A)  - cos (Pi)))<eps,true);
    EXPECT(Max(abs(tan(A)  - tan (Pi)))<eps,true);
    EXPECT(Max(abs(asin(B) -asin (0.25)))<eps,true);
    EXPECT(Max(abs(acos(B) -acos (0.25)))<eps,true);
    EXPECT(Max(abs(atan(A) -atan (Pi)))<eps,true);
    EXPECT(Max(abs(sinh(A) -sinh (Pi)))<eps,true);
    EXPECT(Max(abs(cosh(A) -cosh (Pi)))<eps,true);
    EXPECT(Max(abs(tanh(A) -tanh (Pi)))<eps,true);
    EXPECT(Max(abs(log(A)  -log  (Pi)))<eps,true);
    EXPECT(Max(abs(log10(A)-log10(Pi)))<eps,true);
    EXPECT(Max(abs(exp(A)  -exp  (Pi)))<eps,true);
    EXPECT(Max(abs(sqrt(A) -sqrt (Pi)))<eps,true);
    EXPECT(Max(abs(pow2(A) -    Pi*Pi))<eps,true);
    EXPECT(Max(abs(pow3(A) - Pi*Pi*Pi))<eps,true);
    EXPECT(Max(abs(pow4(A)  -Pi*Pi*Pi*Pi))<100*eps,true);
    //EXPECT(Max(abs(pow10(A)  -pow10(Pi)))<1000*eps,true);
    EXPECT(Max(abs(pow(A,0.6)-pow   (Pi,0.6)))<eps,true);
    EXPECT(Max(abs(pow(A,B)  -pow  (Pi,0.25)))<eps,true);
    EXPECT(Max(abs(atan2(A,B)-atan2(Pi,0.25)))<eps,true);
  }
//##########################################################################
//
//  Test the sorting routines
//
  {
    List<double> A,B;
    for (int i=0;i<1000;i++) A.push_back(i);
    for (int i=0;i<1000;i++) B.push_back(i);
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

  return pass ? 0 : -1;

}
