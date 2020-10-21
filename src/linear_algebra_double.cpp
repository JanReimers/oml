#include "oml/matrix.h"
//---------------------------------------------------------------------------------
//
//  Make template instance
//
#define TYPE  double
#define MTYPE Matrix
#include "src/cnumeric/epsilon.cpp"
//#include "src/numeric/cholsky.cpp"
#include "src/numeric/eigsort.cpp"
//#include "src/numeric/gaussj.cpp"
//#include "src/numeric/invmat.cpp"
//#include "src/numeric/invtri.cpp"
//#include "src/numeric/lubksb.cpp"
//#include "src/numeric/lubksbm.cpp"
//#include "src/numeric/ludcmp.cpp"
//#include "src/numeric/lusolver.cpp"
#include "src/numeric/qldecomp.cpp"
//#include "src/numeric/solvetri.cpp"
#include "src/numeric/svbksb.cpp"
#include "src/numeric/svdcmp.cpp"
#include "src/numeric/svsolver.cpp"
#include "src/numeric/tridiag.cpp"



