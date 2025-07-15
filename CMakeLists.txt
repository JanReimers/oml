cmake_minimum_required(VERSION 3.30.0)
project(oml VERSION 1.1.0 LANGUAGES CXX)

enable_language (Fortran)
set(CXX_SCAN_FOR_MODULES on)
 
add_library(oml "") 
add_library(omlm "") 
set_property(TARGET oml PROPERTY CXX_STANDARD 20)
set_property(TARGET omlm PROPERTY CXX_STANDARD 20)
add_subdirectory(src)
# target_sources(oml PUBLIC )

# -Wno-unknown-pragmas
target_compile_options(oml PRIVATE -fopenmp -Wall -Wno-unknown-pragmas -DHAVE_STDINT_H
    $<$<CONFIG:Debug>:-g -O0>
    $<$<CONFIG:Release>: -O2>
)
target_compile_definitions(oml PRIVATE 
    $<$<CONFIG:Debug>:-DDEBUG>
    $<$<CONFIG:Release>:-DNDEBUG>
)
target_include_directories(oml PRIVATE ./ include )

target_compile_options(omlm PRIVATE -Wall -Wno-unknown-pragmas -DHAVE_STDINT_H
    $<$<CONFIG:Debug>:-g -O0>
    $<$<CONFIG:Release>: -O2>
)
target_compile_definitions(omlm PRIVATE 
    $<$<CONFIG:Debug>:-DDEBUG>
    $<$<CONFIG:Release>:-DNDEBUG>
)
target_include_directories(omlm PRIVATE include )

add_subdirectory(UnitTests)