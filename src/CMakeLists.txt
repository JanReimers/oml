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


 add_subdirectory(linear_algebra)
