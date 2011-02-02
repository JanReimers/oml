#include <cmath>

// Modifications for C++ and oml containers Copyright (1994-2003), Jan N. Reimers

//-------------------------------------------------------------
//
//     estimate unit roundoff in quantities of size x.
//
//    this program should function properly on all systems
//    satisfying the following two assumptions,
//       1.  the base used in representing floating point
//           numbers is not a power of three.
//       2.  the quantity  a  in statement 10 is represented to 
//           the accuracy used in floating point variables
//           that are stored in memory.
//    the statement number 10 and the go to 10 are intended to
//    force optimizing compilers to generate code satisfying 
//    assumption 2.
//    under these assumptions, it should be true that,
//           a  is not exactly equal to four-thirds,
//           b  has a zero for its last bit or digit,
//           c  is not exactly equal to one,
//           eps  measures the separation of 1.0 from
//                the next larger floating point number.
//
double epsilon (double x)
{
  double a,b,c,eps;
  a = 4.0/3.0;
  do
    {
      b = a - 1.0;
      c = b + b + b;
      eps = fabs(c-1.0);
    } while (eps == 0.0);
	
  return eps*fabs(x);
}
