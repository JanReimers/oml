#include "gtest/gtest.h"

int TestEigenDouble();

int main(int argc, char **argv)
{
    testing::InitGoogleTest(&argc, argv);
//    testing::GTEST_FLAG(filter) = "My/MatrixTests/*.BinaryOps";
    return RUN_ALL_TESTS();
}


