# add_executable(UToml main.cc Benchmarks.cpp LinearAlgebraTests.C )

add_executable(UTomlm main.cc UTVector.cpp UTVector3D.cpp  UTMatrix.cpp UTSMatrix.cpp UTDiagonalMatrix.cpp)
set_property(TARGET UTomlm PROPERTY CXX_STANDARD 20)

target_link_options(UTomlm PRIVATE )
target_include_directories(UTomlm PRIVATE ../include )
target_link_libraries(UTomlm omlm gtest)