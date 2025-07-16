target_sources(omlm PUBLIC FILE_SET all_my_modules TYPE CXX_MODULES FILES 
    stream.cpp
    fakedouble.cpp;
    veclimitsm.cpp 
    matlimits.cpp
    shape.cpp
    unop.cpp
    mixtypes.cpp
    xpr.cpp
    indexable.cpp
    arrindex.cpp
    vecindex.cpp
    cow.cpp
    random.cpp
    vector3d.cpp
    matrix3d.cpp
    vectorm.cpp 
    matrixbase.cpp
    matindex.cpp
    rowcol.cpp
    matrixm.cpp
    smatindex.cpp
    smatrix.cpp
    diagonalmatrix.cpp
    numeric.cpp
    imp_numeric.cpp
    indexsort.cpp
    )

target_sources(omlm PRIVATE )

