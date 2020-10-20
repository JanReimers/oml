// FIle: fft.cc  Fast Fourier Transform routines.

// This code comes from netlib.

void rffti(int n, float *wsave, int *ifac)
{
  if (n == 1) return;
  drfti1(n, std::wsave+n, ifac);
}
//
//  forward transform of a real periodic sequence.
//
//  n     the length of the array r to be transformed.  the method
//        is most efficient when n is a product of small primes.
//        n may change so long as different work arrays are provided
//
//  r     a real array of length n which contains the sequence
//        to be transformed
//  std::wsave a work array which must be dimensioned at least 2*n+15.
//        in the program that calls rfftf. the std::wsave array must be
//        initialized by calling subroutine rffti(n,wsave) and a
//        different std::wsave array must be used for each different
//        value of n. this initialization does not have to be
//        repeated so long as n remains unchanged thus subsequent
//        transforms can be obtained faster than the first.
//        the same std::wsave array can be used by rfftf and rfftb.
//
//
//  output parameters
//
//  r     r(1) = the sum from i=1 to i=n of r(i)
//
//        if n is even set l =n/2   , if n is odd set l = (n+1)/2
//
//          then for k = 2,...,l
//
//             r(2*k-2) = the sum from i = 1 to i = n of
//
//                  r(i)*cos((k-1)*(i-1)*2*pi/n)
//
//             r(2*k-1) = the sum from i = 1 to i = n of
//
//               -r(i)*sin((k-1)*(i-1)*2*pi/n)
//
//        if n is even
//
//             r(n) = the sum from i = 1 to i = n of
//
//                  (-1)**(i-1)*r(i)
//
// *****  note
//             this transform is unnormalized since a call of rfftf
//             followed by a call of rfftb will multiply the input
//             sequence by n.
//
//  std::wsave   contains results which must not be destroyed between
//          calls of rfftf or rfftb.
//
void fdrfftf(int n,float *r,float *wsave,int *ifac)
{
  if(n==1)return;
  drftf1(n,r,wsave,wsave+n,ifac);
}

