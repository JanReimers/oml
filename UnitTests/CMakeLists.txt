

add_executable(UTMain main.cc Benchmarks.cpp DiagonalMatrixTests.cpp Epsilons.C LinearAlgebraTests.C MatrixTests.cc SymmetricMatrixTests.cpp UT-Vector3D.cc VectorTests.cc)
target_compile_options(UTMain PRIVATE -fopenmp -Wall -Wno-unknown-pragmas
    $<$<CONFIG:Debug>:-g -O0>
    $<$<CONFIG:Release>: -O2>
)

target_link_options(UTMain PRIVATE -fopenmp)
target_include_directories(UTMain PRIVATE ../include )
target_link_libraries(UTMain oml lapack gtest)