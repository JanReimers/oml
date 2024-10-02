#include "gtest/gtest.h"

int main(int argc, char **argv)
{
    testing::InitGoogleTest(&argc, argv);


//    testing::GTEST_FLAG(filter) = "My/MatrixTests/*.Fill_SetLimits_SubMatrix";
//    testing::GTEST_FLAG(filter) = "MatrixComplexTests.RangeBasedLoops";
//    testing::GTEST_FLAG(filter) = "My/BenchmarkTests/*.*:BenchmarkRealTests.*";
//    testing::GTEST_FLAG(filter) = "BenchmarkRealTests.*";
//    testing::GTEST_FLAG(filter) = "My/DiagonalMatrixTests/*.MatrixAlgebra";
//    testing::GTEST_FLAG(filter) = "DiagonalMatrixComplexTests.MixedTypes";
//    testing::GTEST_FLAG(filter) = "My/SMatrixTests/*.Constructors";
//    testing::GTEST_FLAG(filter) = "My/SMatrixTests/*.Fill_SetLimits_SubMatrix";
//    testing::GTEST_FLAG(filter) = "SMatrixDoubleTests.*";
//    testing::GTEST_FLAG(filter) = "Vector3DTests.*";
//    testing::GTEST_FLAG(filter) = "My/MatrixTests/*.OverloadedOperators1";
    testing::GTEST_FLAG(filter) = "LinearAlgebraTests.*";
    return RUN_ALL_TESTS();
}


