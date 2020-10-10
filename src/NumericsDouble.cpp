#include "oml/matrix.h"
//---------------------------------------------------------------------------------
//
//  Make template instance
//
#define TYPE  double
#define MTYPE Matrix
#include "oml/numeric/cholsky.cc"
#include "oml/numeric/eigsort.cc"
#include "oml/numeric/gaussj.cc"
#include "oml/numeric/invmat.cc"
#include "oml/numeric/invtri.cc"
#include "oml/numeric/lubksb.cc"
#include "oml/numeric/lubksbm.cc"
#include "oml/numeric/ludcmp.cc"
#include "oml/numeric/lusolver.cc"
#include "oml/numeric/qldecomp.cc"
#include "oml/numeric/solvetri.cc"
#include "oml/numeric/svbksb.cc"
#include "oml/numeric/svdcmp.cc"
#include "oml/numeric/svsolver.cc"
#include "oml/numeric/tridiag.cc"



