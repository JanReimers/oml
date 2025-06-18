

add_executable(UToml main.cc Benchmarks.cpp DiagonalMatrixTests.cpp Epsilons.C LinearAlgebraTests.C MatrixTests.cc SymmetricMatrixTests.cpp UT-Vector3D.cc VectorTests.cc)
target_compile_options(UToml PRIVATE -fopenmp -Wall -Wno-unknown-pragmas
    $<$<CONFIG:Debug>:-g -O0>
    $<$<CONFIG:Release>: -O2>
)

target_link_options(UToml PRIVATE -fopenmp)
target_include_directories(UToml PRIVATE ../include )
target_link_libraries(UToml oml lapack gtest)