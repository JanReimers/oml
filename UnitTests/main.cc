#include "gtest/gtest.h"

int TestEigenDouble();

int main(int argc, char **argv)
{
    TestEigenDouble();
    testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}


