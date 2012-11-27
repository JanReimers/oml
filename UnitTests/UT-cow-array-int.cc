// File: UT-cow-array-int.cc  Unit test the Array class for int data types.

// Copyright (1994-2003), Jan N. Reimers

#include "gtest/gtest.h"
#include "oml/UnitTests/UnitTest.H"
#include "oml/cow.h"

//
//  Simple tests to allow valgrind to memory leaks
//
class COWArrayTesting : public ::testing::Test
{
public:
    COWArrayTesting() {};
};

TEST_F(COWArrayTesting,SizeConstructor)
{
    cow_array<int> A1(10);
}

TEST_F(COWArrayTesting,ShallowCopyonstructor)
{
    cow_array<int> A1(10);
    cow_array<int> A2(A1);
}

TEST_F(COWArrayTesting,ShallowCopyConstructor1)
{
    cow_array<int> A1(10);
    {
        cow_array<int> A2(A1);
    }
}


TEST_F(COWArrayTesting,ShallowCopyConstructor2)
{
    cow_array<int>* A1 =new cow_array<int> (10);
    cow_array<int> A2(*A1);
    delete A1;
}

TEST_F(COWArrayTesting,DeepCopyonstructor)
{
    cow_array<int> A1(10);
    cow_array<int> A2(A1);
    A1.Get()[5]=4;
}

TEST_F(COWArrayTesting,DeepCopyConstructor1)
{
    cow_array<int> A1(10);
    {
        cow_array<int> A2(A1);
    }
    A1.Get()[5]=4;
}


TEST_F(COWArrayTesting,DeepCopyConstructor2)
{
    cow_array<int>* A1 =new cow_array<int> (10);
    cow_array<int> A2(*A1);
    A1->Get()[5]=4;
    delete A1;
}

TEST_F(COWArrayTesting,ShallowOpEq)
{
    cow_array<int> A1(10);
    cow_array<int> A2(10);
    A2=A1;
}

TEST_F(COWArrayTesting,ShallowCopyOpEq1)
{
    cow_array<int> A1(10);
    {
        cow_array<int> A2(10);
        A2=A1;
    }
}
TEST_F(COWArrayTesting,ShallowCopyOpEq2)
{
    cow_array<int> A1(10);
    cow_array<int> A2(10);
    A1=A2;
}



TEST_F(COWArrayTesting,ShallowCopyOpEq3)
{
    cow_array<int>* A1 =new cow_array<int> (10);
    cow_array<int> A2(10);
    A2=*A1;
    delete A1;
}

TEST_F(COWArrayTesting,ShallowCopyOpEq4)
{
    cow_array<int>* A1 =new cow_array<int> (10);
    cow_array<int> A2(10);
    *A1=A2;
    delete A1;
}


TEST_F(COWArrayTesting,DeepCopyOpEq1)
{
    cow_array<int> A1(10);
    {
        cow_array<int> A2(10);
        A2=A1;
        A1.Get()[5]=4;
    }
}
TEST_F(COWArrayTesting,DeepCopyOpEq2)
{
    cow_array<int> A1(10);
    cow_array<int> A2(10);
    A1=A2;
    A1.Get()[5]=4;

}



TEST_F(COWArrayTesting,DeepCopyOpEq3)
{
    cow_array<int>* A1 =new cow_array<int> (10);
    cow_array<int> A2(10);
    A2=*A1;
    A2.Get()[5]=4;

    delete A1;
}

TEST_F(COWArrayTesting,DeepCopyOpEq4)
{
    cow_array<int>* A1 =new cow_array<int> (10);
    cow_array<int> A2(10);
    *A1=A2;
    A2.Get()[5]=4;
    delete A1;
}
