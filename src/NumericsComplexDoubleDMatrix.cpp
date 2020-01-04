#include "oml/dmatrix.h"
#include <complex>
//---------------------------------------------------------------------------------
//
//  Make template instance
//
#define TYPE  std::complex<double>
#define MTYPE DMatrix
#include "numeric/svdcmp.cc"



