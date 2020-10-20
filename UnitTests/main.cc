#include "gtest/gtest.h"

int TestEigenDouble();

int main(int argc, char **argv)
{
    testing::InitGoogleTest(&argc, argv);
//    testing::GTEST_FLAG(filter) = "My/MatrixTests/*.BinaryOps";
//    testing::GTEST_FLAG(filter) = "MatrixComplexTests.RangeBasedLoops";
//    testing::GTEST_FLAG(filter) = "MatrixRealTests.MinMax";
    return RUN_ALL_TESTS();
}


