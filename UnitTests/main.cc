#include "gtest/gtest.h"

int TestEigenDouble();

int main(int argc, char **argv)
{
    testing::InitGoogleTest(&argc, argv);


//    testing::GTEST_FLAG(filter) = "My/MatrixTests/*.AsciiAndBinaryIO";
//    testing::GTEST_FLAG(filter) = "My/MatrixTests/*.BinaryOps";
//    testing::GTEST_FLAG(filter) = "MatrixComplexTests.RangeBasedLoops";
//    testing::GTEST_FLAG(filter) = "My/DiagonalMatrixTests/*.MatrixAlgebra";
//    testing::GTEST_FLAG(filter) = "DiagonalMatrixComplexTests.*";
    return RUN_ALL_TESTS();
}


