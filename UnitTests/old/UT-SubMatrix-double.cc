// File: UT-Array-int.cc  Unit test the Array class for int data types.

// Copyright (1994-2003), Jan N. Reimers

#include "gtest/gtest.h"
#include "oml/UnitTests/UnitTest.H"
#include "oml/submatrix.h"
#include <stdlib.h>
#include <fstream>
#include <iostream>

using std::cout;
using std::endl;

class SubMatrixTesting : public ::testing::Test
{
public:
    SubMatrixTesting() :A1(10,10), A2(100,100)
    {
        StreamableObject::SetToPretty();
        MatLimits lim=A1.GetLimits();
        for (int i=lim.Row.Low;i<=lim.Row.High;i++)
        for (int j=lim.Col.Low;j<=lim.Col.High;j++)
        A1(i,j)=100*j+i;
    }

    template <class T,class A,Store M, Data D> inline
    void TestMatrixSubscripting(const Indexable<T,A,M,D,MatrixShape>& m)
    {
        double t=m(5,2); //Matrix subscripting is OK
    }

    Matrix<double> A1,A2;
};

TEST_F(SubMatrixTesting,DefaultConstructor)
{
    cout << "A1=" << A1 << endl;
    RefSubMatrix<double> B1(A1,2,8,3,7);
    TestMatrixSubscripting(B1);
    cout << "B1=" << B1 << endl;
}

