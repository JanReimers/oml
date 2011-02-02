// File: UT-List-int.cc  Unit test the List class for int data types.

// Copyright (1994-2003), Jan N. Reimers

#include "oml/UnitTests/UnitTest.H"
#include "oml/list.h"
#include "oml/list_io.h"
#include "oml/array.h"
#include "oml/array_io.h"
#include "oml/random.h"
#include "oml/minmax.h"
#include <fstream>

typedef List<int> ListType;

int main()
{
  const char* Class="List<int>";
  bool pass=true;
//
//  Start message.
//

  StartClass(Class);
  StreamableObject::SetOutputMode(StreamableObject::pretty);

//###########################################################################
//
//  Constructors
//
  {
    EXPECT(ListType(),"{ }");
    ListType A1,A2;
    for (index_t i=0;i<10;i++) A1.push_back(i);
    EXPECT(A1,"{ 0 1 2 3 4 5 6 7 8 9 }");
    EXPECT(ListType(A1),"{ 0 1 2 3 4 5 6 7 8 9 }");
    EXPECT(A2=A1,"{ 0 1 2 3 4 5 6 7 8 9 }");
    EXPECT2(A1,SetSize(13,true),"{ 0 1 2 3 4 5 6 7 8 9 * * * }");
    EXPECT2(A1,SetSize(10),"{ * * * * * * * * * * }");
    EXPECT3(A2,FillLinear(A2,1,10),"{ 1 2 3 4 5 6 7 8 9 10 }");
    EXPECT(A2[8],"9");
    EXPECT(A2.SubList(4,7),"{ 5 6 7 8 }");
    EXPECT(A2&A2,"{ 1 2 3 4 5 6 7 8 9 10 1 2 3 4 5 6 7 8 9 10 }");
    EXPECT3(A1,A1.Empty(),"{ }");
  }
  {
    ListType A,B,C;
    for (index_t i=0;i<10;i++) A.push_back(i);
    C=B=A;
    C[0]=-1;
    EXPECT(A,"{ 0 1 2 3 4 5 6 7 8 9 }");
    EXPECT(B,"{ 0 1 2 3 4 5 6 7 8 9 }");
    EXPECT(C,"{ -1 1 2 3 4 5 6 7 8 9 }");  
    EXPECT3(A,A.push_back(999),"{ 0 1 2 3 4 5 6 7 8 9 999 }");
    EXPECT3(A,A.InsertAt(-888,1),"{ 0 -888 1 2 3 4 5 6 7 8 9 999 }");
    EXPECT(A.find(4),5);
    EXPECT(A.find(10),A.size());    
    EXPECT3(A,A.Remove(5),"{ 0 -888 1 2 3 4 6 7 8 9 999 }");
    A.push_back(6);
    A.InsertAt(6,5);
    EXPECT3(A,A.RemoveAll(6),"{ 0 -888 1 2 3 4 7 8 9 999 }");
    EXPECT3(A,A.RemoveAt(0),"{ -888 1 2 3 4 7 8 9 999 }");
    EXPECT(A.front(),-888);
    EXPECT(A.back(),999); 
    Array<index_t> index(A.size());
    FillLinear(index,A.size()-1,0);
    EXPECT3(A,A.ReIndex(index),"{ 999 9 8 7 4 3 2 1 -888 }");
  }
  //###########################################################################
//
//  Algebra overloaded operator stuff.
//
  {
    ListType A,B;
    for (int i=0;i<6;i++) 
    {
      A.push_back(13);
      B.push_back(5);
    }
    EXPECT(A+B,"{ 18 18 18 18 18 18 }");
    EXPECT(A-B,"{ 8 8 8 8 8 8 }");
    EXPECT(A*B,"{ 65 65 65 65 65 65 }");
    EXPECT(A/B,"{ 2 2 2 2 2 2 }");            
    EXPECT(A+2,"{ 15 15 15 15 15 15 }");
    EXPECT(A-2,"{ 11 11 11 11 11 11 }");
    EXPECT(A*2,"{ 26 26 26 26 26 26 }");
    EXPECT(A/2,"{ 6 6 6 6 6 6 }");
    EXPECT(A+=B,"{ 18 18 18 18 18 18 }");
    EXPECT(A-=B,"{ 13 13 13 13 13 13 }");
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
    //  Verify that mixing type also works.
    //
    EXPECT(A+2.1,"{ 15.1 15.1 15.1 15.1 15.1 15.1 }");
    EXPECT(2.1+A,"{ 15.1 15.1 15.1 15.1 15.1 15.1 }");
    EXPECT(A-2.1,"{ 10.9 10.9 10.9 10.9 10.9 10.9 }");
    EXPECT(2.1-A,"{ -10.9 -10.9 -10.9 -10.9 -10.9 -10.9 }");

    EXPECT3(A,FillRandom(A,100),"{ * * * * * * }");
    EXPECT(Max(A)<=100,"1");
    EXPECT(Max(A)>=0,"1");
    
    {
      ListType Ran;
      for (int i=0;i<1000;i++) Ran.push_back(i);
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
  return pass ? 0 : -1;
}



    
