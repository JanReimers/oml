

add_executable(UToml main.cc Benchmarks.cpp DiagonalMatrixTests.cpp Epsilons.C LinearAlgebraTests.C MatrixTests.cc SymmetricMatrixTests.cpp UT-Vector3D.cc VectorTests.cc)
target_compile_options(UToml PRIVATE -Wall -Wno-unknown-pragmas
    $<$<CONFIG:Debug>:-g -O0>
    $<$<CONFIG:Release>: -O2>
)

target_link_options(UToml PRIVATE -fopenmp)
target_include_directories(UToml PRIVATE ../include )
target_link_libraries(UToml oml lapack blas gtest)

add_executable(UTomlm main.cc UTvectorm.cpp UTvector3Dm.cpp  UTmatrixm.cpp)
set_property(TARGET UTomlm PROPERTY CXX_STANDARD 20)

target_link_options(UTomlm PRIVATE )
target_include_directories(UTomlm PRIVATE ../include )
target_link_libraries(UTomlm omlm gtest)