// File: UT-Array-int.cc  Unit test the Array class for int data types.

// Copyright (1994-2003), Jan N. Reimers

#include "oml/UnitTests/UnitTest.H"
#include "oml/array.h"
#include "oml/array_io.h"
#include "oml/minmax.h"
#include "oml/random.h"
#include <stdlib.h>
#include <fstream>


int main()
{
  const char* Class="Array<int>";
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
      EXPECT(Array<index_t>(),"{ }");
      EXPECT(Array<index_t>(10),"{ * * * * * * * * * * }");
      
      Array<index_t> A1(10),A2(10);
      Fill(A1,3);
      EXPECT(Array<index_t>(A1),"{ 3 3 3 3 3 3 3 3 3 3 }");      
      EXPECT(A2=A1,"{ 3 3 3 3 3 3 3 3 3 3 }");      
      EXPECT2(A1,SetSize(13,true),"{ 3 3 3 3 3 3 3 3 3 3 * * * }");
      EXPECT2(A1,SetSize(10),"{ * * * * * * * * * * }");
      EXPECT3(A2,FillLinear(A2,1,10),"{ 1 2 3 4 5 6 7 8 9 10 }");
      EXPECT(A2[8],"9");
      EXPECT(A2.SubArray(4,7),"{ 5 6 7 8 }");
      EXPECT(A2&A2,"{ 1 2 3 4 5 6 7 8 9 10 1 2 3 4 5 6 7 8 9 10 }");
   }

//###########################################################################
//
//  Algebra overloaded operator stuff.
//
  {
    Array<index_t> A(6),B(6);
    Fill(A,13);
    Fill(B, 5);
    EXPECT(A+B,"{ 18 18 18 18 18 18 }");
    EXPECT(A-B,"{ 8 8 8 8 8 8 }");
    EXPECT(A*B,"{ 65 65 65 65 65 65 }");
    EXPECT(A/B,"{ 2 2 2 2 2 2 }");            
    EXPECT(A+2,"{ 15 15 15 15 15 15 }");
    EXPECT(A-2,"{ 11 11 11 11 11 11 }");
    EXPECT(A*2,"{ 26 26 26 26 26 26 }");
    EXPECT(A/2,"{ 6 6 6 6 6 6 }");
    EXPECT(2+A,"{ 15 15 15 15 15 15 }");
    EXPECT(2-A,"{ -11 -11 -11 -11 -11 -11 }");
    EXPECT(2*A,"{ 26 26 26 26 26 26 }");
    EXPECT(26/A,"{ 2 2 2 2 2 2 }");
    EXPECT(A+=B,"{ 18 18 18 18 18 18 }");
    EXPECT(A-=B,"{ 13 13 13 13 13 13 }");
    EXPECT(A*=B,"{ 65 65 65 65 65 65 }");
    EXPECT(A/=B,"{ 13 13 13 13 13 13 }");
		EXPECT(A+=2,"{ 15 15 15 15 15 15 }");            
    EXPECT(A-=2,"{ 13 13 13 13 13 13 }");      
    EXPECT(A*=2,"{ 26 26 26 26 26 26 }");      
    EXPECT(A/=2,"{ 13 13 13 13 13 13 }");      
    EXPECT(Dot(A,B),13*5*6);
    EXPECT(A==A,"1");
    EXPECT(A==B,"0");      
    EXPECT(A!=A,"0");
    EXPECT(A!=B,"1");
    EXPECT(Sum(A),"78");
    EXPECT(Sum(pow2(A)),"1014");
    //
    //  Verify that mixing types also works.
    //
    EXPECT(A+2.1,"{ 15.1 15.1 15.1 15.1 15.1 15.1 }");
    EXPECT(2.1+A,"{ 15.1 15.1 15.1 15.1 15.1 15.1 }");
    EXPECT(A-2.1,"{ 10.9 10.9 10.9 10.9 10.9 10.9 }");
    EXPECT(2.1-A,"{ -10.9 -10.9 -10.9 -10.9 -10.9 -10.9 }");

    EXPECT3(A,FillRandom(A,100),"{ * * * * * * }");
    EXPECT(Max(A)<=100,"1");
    EXPECT(Max(A)>=0,"1");
    
    {
      Array<int> Ran(1000);
      FillRandom(Ran,100);
      index_t min=Ran[0], max=Ran[0];
      for (int i=0;i<1000;i++)
      {
	min=Min(min,Ran[i]);
	max=Max(max,Ran[i]);
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
    Array<int> A(1000),B;
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
    Array<int> A(1000),B;
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
