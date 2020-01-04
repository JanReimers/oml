// File: htribk.cc  Back transform real symmetric matrix.
#include "oml/matrix.h"
#include <cassert>
#include <complex>

#ifndef TYPE
#error "htribk.cc TYPE not defined"
#endif

#ifndef MTYPE
#error "htribk.cc MTYPE not defined"
#endif
using std::real;

// Modifications for oml containers Copyright (1994-2003), Jan N. Reimers

//--------------------------------------------------------------------
//
//     this subroutine is a translation of a complex analogue of
//     the algol procedure trbak1, num. math. 11, 181-195(1968)
//     by martin, reinsch, and wilkinson.
//     handbook for auto. comp., vol.ii-linear algebra, 212-226(1971).
//
//    this subroutine forms the eigenvectors of a complex hermitian
//    matrix by back transforming those of the corresponding
//    real symmetric tridiagonal matrix determined by  htridi.
//
template <class T, class M> void htribk(const M& A,const Vector<std::complex<T> >& tau, M& Z)
{
  int i,j,k,l,m,n;
  T h;

  m=A.GetNumRows();
  n=m;
  assert(A.GetNumRows()==A.GetNumCols());
  assert(m==Z.GetNumRows());
  assert(m==Z.GetNumCols());
  assert(m==tau.size());

  typename M::Subscriptor sZ(Z);

  if (m == 0) return;
  for (k=1;k<=n;k++)
    for (j=1;j<=m;j++)
      sZ(k,j) = real(sZ(k,j))*conj(tau(k));

  if (n == 1) return;
  for(i=2;i<=n;i++)
  {
    l = i - 1;
    h = imag(A(i,i));
    if (h == 0.0) break;
    for (j=1;j<=m;j++)
    {
      std::complex<T> s(0.0);
      for (k=1;k<=l;k++) s += A(i,k)*sZ(k,j);
      s = (s / h) / h; // double divisions avoid possible underflow.
      for (k=1;k<=l;k++) sZ(k,j) -= s*conj(A(i,k));
    }
  }
  return;
}

typedef TYPE Type;
typedef MTYPE<std::complex<Type> > MType;
template void htribk<Type,MType>(const MType& A,const Vector<std::complex<Type> >& tau, MType& Z);

