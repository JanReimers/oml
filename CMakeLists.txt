cmake_minimum_required(VERSION 3.20.0)
project(oml VERSION 1.0.0 LANGUAGES CXX)

enable_language (Fortran)
#file (GLOB_RECURSE oml_SOURCES CONFIGURE_DEPENDS "src/*.cpp")
# file (GLOB_RECURSE oml_HEADERS CONFIGURE_DEPENDS "include/*.h")
 
add_library(oml "") 
add_subdirectory(src)
# target_sources(oml PUBLIC )

# -Wno-unknown-pragmas
target_compile_options(oml PRIVATE -Wall -Wno-unknown-pragmas -DHAVE_STDINT_H
    $<$<CONFIG:Debug>:-g -O0>
    $<$<CONFIG:Release>: -O2>
)
target_compile_definitions(oml PRIVATE 
    $<$<CONFIG:Debug>:-DDEBUG>
    $<$<CONFIG:Release>:-DNDEBUG>
)
target_include_directories(oml PRIVATE ./ include )

add_subdirectory(UnitTests)