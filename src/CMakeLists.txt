#cmake_minimum_required(VERSION 3.5.0)
#project(oml VERSION 1.0.0 LANGUAGES CXX)

target_sources(oml PRIVATE 
    stream.cpp
    veclimit.cpp 
    matlimit.cpp                
    ran250.cpp                        
    instance_int.cpp       
    instance_double.cpp      
    instance_matrix_double.cpp  
    instance_complex.cpp  
    InstanceComplexVec3.cpp  
    InstanceVec3.cpp            
    linear_algebra_complex.cpp  
    linear_algebra_double.cpp   
 )

 #target_sources(omlm PRIVATE )
 target_sources(omlm PUBLIC FILE_SET all_my_modules TYPE CXX_MODULES FILES vectorm.cpp matrixm.cpp)


 add_subdirectory(linear_algebra)
