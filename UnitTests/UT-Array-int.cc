// File: UT-Array-int.cc  Unit test the Array class for int data types.

// Copyright (1994-2003), Jan N. Reimers

#include "gtest/gtest.h"
#include "oml/UnitTests/UnitTest.H"
#include "oml/array.h"
#include "oml/array_io.h"
#include "oml/minmax.h"
#include "oml/random.h"
#include <stdlib.h>
#include <fstream>

class ArrayTesting : public ::testing::Test
{
public:
    ArrayTesting() :A1(10), A2(10)
    {
        StreamableObject::SetToPretty();
        Fill(A1,3);
    }

    Array<index_t> A1,A2;
};

TEST_F(ArrayTesting,DefaultConstructor)
{
    EXPECT_EQ(ToString(Array<index_t>()),"{ }");
}


TEST_F(ArrayTesting,SizeConstructor)
{
    bool b=EXPECT(Array<index_t>(10),"{ * * * * * * * * * * }");
    EXPECT_TRUE(b);
}

TEST_F(ArrayTesting,CopyConstructor)
{
    EXPECT_EQ(ToString(Array<index_t>(A1)),"{ 3 3 3 3 3 3 3 3 3 3 }");
}

TEST_F(ArrayTesting,OpEq)
{
    EXPECT_EQ(ToString(A2=A1),"{ 3 3 3 3 3 3 3 3 3 3 }");
}

TEST_F(ArrayTesting,DefaulConstructor)
{
    EXPECT2(A1,SetSize(13,true)," { 3 3 3 3 3 3 3 3 3 3 * * * }");
    EXPECT2(A1,SetSize(10,true)," { * * * * * * * * * * }");
}

TEST_F(ArrayTesting,FillSubscriptSubArrayAndConcatenation)
{
    EXPECT3(A2,FillLinear(A2,1,10)," { 1 2 3 4 5 6 7 8 9 10 }");
    EXPECT_EQ(ToString(A2[8]),"9");
    EXPECT_EQ(ToString(A2.SubArray(4,7)),"{ 5 6 7 8 }");
    EXPECT_EQ(ToString(A2&A2),"{ 1 2 3 4 5 6 7 8 9 10 1 2 3 4 5 6 7 8 9 10 }");
}

TEST_F(ArrayTesting,OverloadedOperators)
{
    Array<index_t> A(6),B(6);
    Fill(A,13);
    Fill(B, 5);
    EXPECT_TRUE(EXPECT(A+B," { 18 18 18 18 18 18 }"));
    EXPECT_TRUE(EXPECT(A-B," { 8 8 8 8 8 8 }"));
    EXPECT_TRUE(EXPECT(A*B," { 65 65 65 65 65 65 }"));
    EXPECT_TRUE(EXPECT(A/B," { 2 2 2 2 2 2 }"));
    EXPECT_TRUE(EXPECT(A+2," { 15 15 15 15 15 15 }"));
    EXPECT_TRUE(EXPECT(A-2," { 11 11 11 11 11 11 }"));
    EXPECT_TRUE(EXPECT(A*2," { 26 26 26 26 26 26 }"));
    EXPECT_TRUE(EXPECT(A/2," { 6 6 6 6 6 6 }"));
    EXPECT_TRUE(EXPECT(2+A," { 15 15 15 15 15 15 }"));
    EXPECT_TRUE(EXPECT(2-A," { -11 -11 -11 -11 -11 -11 }"));
    EXPECT_TRUE(EXPECT(2*A," { 26 26 26 26 26 26 }"));
    EXPECT_TRUE(EXPECT(26/A," { 2 2 2 2 2 2 }"));
    EXPECT_TRUE(EXPECT(A+=B," { 18 18 18 18 18 18 }"));
    EXPECT_TRUE(EXPECT(A-=B," { 13 13 13 13 13 13 }"));
    EXPECT_TRUE(EXPECT(A*=B," { 65 65 65 65 65 65 }"));
    EXPECT_TRUE(EXPECT(A/=B," { 13 13 13 13 13 13 }"));
    EXPECT_TRUE(EXPECT(A+=2," { 15 15 15 15 15 15 }"));
    EXPECT_TRUE(EXPECT(A-=2," { 13 13 13 13 13 13 }"));
    EXPECT_TRUE(EXPECT(A*=2," { 26 26 26 26 26 26 }"));
    EXPECT_TRUE(EXPECT(A/=2," { 13 13 13 13 13 13 }"));
    EXPECT_TRUE(EXPECT(Dot(A,B),13*5*6));
    EXPECT_TRUE(EXPECT(A==A,"1"));
    EXPECT_TRUE(EXPECT(A==B,"0"));
    EXPECT_TRUE(EXPECT(A!=A,"0"));
    EXPECT_TRUE(EXPECT(A!=B,"1"));
    EXPECT_TRUE(EXPECT(Sum(A),"78"));
    EXPECT_TRUE(EXPECT(Sum(pow2(A)),"1014"));
    //
    //  Verify that mixing types also works.
    //
    EXPECT_TRUE(EXPECT(A+2.1," { 15.1 15.1 15.1 15.1 15.1 15.1 }"));
    EXPECT_TRUE(EXPECT(2.1+A," { 15.1 15.1 15.1 15.1 15.1 15.1 }"));
    EXPECT_TRUE(EXPECT(A-2.1," { 10.9 10.9 10.9 10.9 10.9 10.9 }"));
    EXPECT_TRUE(EXPECT(2.1-A," { -10.9 -10.9 -10.9 -10.9 -10.9 -10.9 }"));

    EXPECT3(A,FillRandom(A,100)," { * * * * * * }");
    EXPECT_TRUE(EXPECT(Max(A)<=100,"1"));
    EXPECT_TRUE(EXPECT(Max(A)>=0,"1"));

    {
        Array<int> Ran(1000);
        FillRandom(Ran,100);
        index_t min=Ran[0], max=Ran[0];
        for (int i=0; i<1000; i++)
        {
            min=Min(min,Ran[i]);
            max=Max(max,Ran[i]);
        }
        EXPECT_TRUE(EXPECT(Min(Ran),min));
        EXPECT_TRUE(EXPECT(Max(Ran),max));
    }
}


TEST_F(ArrayTesting,AsciiIO)
{
    StreamableObject::SetToAscii();
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
    EXPECT_EQ(B,A);
}
TEST_F(ArrayTesting,BinaryIO)
{
    StreamableObject::SetToBinary();
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
    EXPECT_EQ(B,A);
}



