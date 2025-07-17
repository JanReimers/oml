# add_executable(UToml main.cc Benchmarks.cpp LinearAlgebraTests.C )

add_executable(UTomlm main.cc UTVector3D.cpp UTMatrix3D.cpp UTVector.cpp UTMatrix.cpp UTSMatrix.cpp UTDiagonalMatrix.cpp UTLinearAlgebra.C Epsilons.C)
set_property(TARGET UTomlm PROPERTY CXX_STANDARD 20)

target_link_options(UTomlm PRIVATE )
target_include_directories(UTomlm PRIVATE ../include )
target_link_libraries(UTomlm omlm gtest lapack blas)
target_compile_options(UTomlm PRIVATE -Wall -Wno-unknown-pragmas 
    $<$<CONFIG:Debug>:-g -O0>
    $<$<CONFIG:Release>: -O2>
)
target_compile_definitions(UTomlm PRIVATE 
    $<$<CONFIG:Debug>:-DDEBUG>
    $<$<CONFIG:Release>:-DNDEBUG>
)