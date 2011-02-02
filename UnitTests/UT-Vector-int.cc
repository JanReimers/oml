// File: UT-Vector-int.cc  Unit test the Vector class for int data types.

// Copyright (1994-2003), Jan N. Reimers

#include "oml/UnitTests/UnitTest.H"
#include "oml/vector.h"
#include "oml/vector_io.h"
#include "oml/minmax.h"
#include "oml/random.h"
#include <fstream>

int main()
{
  const char* Class="Vector<int>";
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
      std::cout << "Constructor tests" << std::endl;
      EXPECT(Vector<index_t>(),"(1:0){ }");
      EXPECT(Vector<index_t>(10),"(1:10){ * * * * * * * * * * }");
      EXPECT(Vector<index_t>(5,14),"(5:14){ * * * * * * * * * * }");      
      EXPECT(Vector<index_t>(VecLimits(5,14)),"(5:14){ * * * * * * * * * * }");      
      
      Vector<index_t> A1(5,14),A2(5,14);
      Fill(A1,3);
      EXPECT(Vector<index_t>(A1),"(5:14){ 3 3 3 3 3 3 3 3 3 3 }");      
      EXPECT(A2=A1,"(5:14){ 3 3 3 3 3 3 3 3 3 3 }");      
      EXPECT2(A1,SetLimits(VecLimits(-2,17),true),"(-2:17){ * * * * * * * 3 3 3 3 3 3 3 3 3 3 * * * }");
      EXPECT2(A1,SetLimits(VecLimits(-6,13)),"(-6:13){ * * * * * * * * * * * * * * * * * * * * }");
      EXPECT3(A2,FillLinear(A2,1,10),"(5:14){ 1 2 3 4 5 6 7 8 9 10 }");
      EXPECT(A2(8),"4");
      EXPECT(A2.SubVector(9,12),"(9:12){ 5 6 7 8 }");
      EXPECT(A2&A2,"(5:24){ 1 2 3 4 5 6 7 8 9 10 1 2 3 4 5 6 7 8 9 10 }");
   }

//###########################################################################
//
//  std::vector + - * / std::vector.
//
  {
    Vector<index_t> A(-3,2),B(-3,2);
    Fill(A,13);
    Fill(B, 5);
    
    EXPECT(A+B,"(-3:2){ 18 18 18 18 18 18 }");
    EXPECT(A-B,"(-3:2){ 8 8 8 8 8 8 }");
    EXPECT(A*B,"390");
    EXPECT(A+2,"(-3:2){ 15 15 15 15 15 15 }");
    EXPECT(A-2,"(-3:2){ 11 11 11 11 11 11 }");
    EXPECT(A*2,"(-3:2){ 26 26 26 26 26 26 }");
    EXPECT(A/2,"(-3:2){ 6 6 6 6 6 6 }");
    EXPECT(A+=B,"(-3:2){ 18 18 18 18 18 18 }");
    EXPECT(A-=B,"(-3:2){ 13 13 13 13 13 13 }");
    EXPECT(A+=2,"(-3:2){ 15 15 15 15 15 15 }");            
    EXPECT(A-=2,"(-3:2){ 13 13 13 13 13 13 }");      
    EXPECT(A*=2,"(-3:2){ 26 26 26 26 26 26 }");      
    EXPECT(A/=2,"(-3:2){ 13 13 13 13 13 13 }");      
    EXPECT(DirectMultiply(A,B),"(-3:2){ 65 65 65 65 65 65 }");      
    EXPECT(DirectDivide(A,B),"(-3:2){ 2 2 2 2 2 2 }");            
    EXPECT(A==A,"1");
    EXPECT(A==B,"0");      
    EXPECT(A!=A,"0");
    EXPECT(A!=B,"1");
    EXPECT(Sum(A),"78");
    EXPECT(Sum(pow2(A)),"1014");
    //
    //  Verify that mixing type also works.
    //
    EXPECT(A+2.1,"(-3:2){ 15.1 15.1 15.1 15.1 15.1 15.1 }");
    EXPECT(2.1+A,"(-3:2){ 15.1 15.1 15.1 15.1 15.1 15.1 }");
    EXPECT(A-2.1,"(-3:2){ 10.9 10.9 10.9 10.9 10.9 10.9 }");
    EXPECT(2.1-A,"(-3:2){ -10.9 -10.9 -10.9 -10.9 -10.9 -10.9 }");

    EXPECT3(A,FillRandom(A,100),"(-3:2){ * * * * * * }");
    EXPECT(Max(A)<= 100,"1");
    EXPECT(Min(A)>=-100,"1");
    
    {
      Vector<int> Ran(1000);
      FillRandom(Ran,100);
      index_t min=Ran(1), max=Ran(1);
      for (int i=1;i<=1000;i++)
      {
	min=Min(min,Ran(i));
	max=Max(max,Ran(i));
      }
      EXPECT(Min(Ran),min);
      EXPECT(Max(Ran),max);
    }
  }
  
//-----------------------------------------------------
//
//  IO tests
//
  StreamableObject::SetToAscii();
  {
    Vector<int> A(-10,1000),B;
    FillRandom(A,100);
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
    Vector<int> A(-10,1000),B;
    FillRandom(A,100);
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
