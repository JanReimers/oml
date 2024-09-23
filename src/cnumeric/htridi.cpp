// File: htridi.cc  Reduce complex hermitian matrix to real tri diagonal form.
#include "oml/vector.h"
#include "oml/matrix.h"
#include <cassert>
#include <cmath>
#include <complex>

#ifndef TYPE
#error "htridi.cc TYPE not defined"
#endif

#ifndef MTYPE
#error "htridi.cc MTYPE not defined"
#endif

// Modifications for oml containers Copyright (1994-2003), Jan N. Reimers



//----------------------------------------------------------------------
//
//    this subroutine reduces a complex hermitian matrix
//    to a real symmetric tridiagonal matrix using
//    unitary similarity transformations.
//
template <class T, class M> void htridi(M& A, Vector<T>& d ,Vector<T>& evals,
	    Vector<T>& e2,Vector<std::complex<T> >& tau)
{
  assert(d.size()==evals .size());
  assert(d.size()==e2.size());
  assert(d.size()==A.GetNumRows());
  assert(d.size()==A.GetNumCols());
  assert(d.size()==tau.size());

  Vector<std::complex<T> > ee(tau.GetLimits());

  typename M::Subscriptor sA(A);
  typename Vector<T                >::Subscriptor sd(d);
  typename Vector<std::complex<T > >::Subscriptor se(ee);
  typename Vector<std::complex<T > >::Subscriptor stau(tau);
  index_t i,j,k,l,ii,jp1;
  double  f,g,h,hh,scale;
  index_t n=d.size();

  stau(n) = se(n) = std::complex<T>(1.0,0.0);
  d=real(A.GetDiagonal());

  for (ii=1;ii<=n;ii++)//    do 300 ii = 1, n
  {
    i = n + 1 - ii;;
    l = i - 1;
    h = 0.0;
    scale = 0.0;
    if (l < 1) goto L130;
    for (k=1;k<=l;k++) scale = scale + fabs(real(sA(i,k))) + fabs(imag(sA(i,k)));

    if (scale != 0.0) goto L140;
    stau(l) = std::complex<T>(1.0,0.0);
		se(l)= 0.0;
    L130:
    e2(i) = 0.0;
    goto L290;

    L140:
    for (k=1;k<=l;k++)
    {
      sA(i,k) /= scale;
      h += norm(sA(i,k));
    }

    e2(i) = scale * scale * h;
    g = sqrt(h);
    se(i) = std::complex<T>(scale * g,imag(se(i)));
    f = pythag(real(sA(i,l)),imag(sA(i,l)));

    if (f == 0.0) goto L160;

    stau(l) = -conj(sA(i,l)*std::complex<T>(real(stau(i)),imag(se(i))))/f;
    h = h + f * g;
    g = 1.0 + g / f;
    sA(i,l) *= g;
    if (l == 1) goto L270;
    goto L170;
L160:
    stau(l) = std::complex<T>(-real(stau(i)),imag(stau(i)));
    sA(i,l) = std::complex<T>(g,imag(sA(i,l)));
L170:
    f = 0.0;

    for (j=1;j<=l;j++)//     do 240 j = 1, l
    {
      std::complex<T> gg(0.0);
      for (k=1;k<=j;k++) gg+=sA(j,k)*conj(sA(i,k));

      jp1 = j + 1;
      if (l < jp1) goto L220;
      for (k=jp1;k<=l;k++) gg+=conj(sA(k,j)*sA(i,k));
L220:
      se(j) = gg / h;
      f +=  real(se(j)*(sA(i,j)));
    }
    hh = f / (h + h);
    for (j=1;j<=l;j++)
    {
      se(j) = conj(se(j)-hh*conj(sA(i,j)));
      for (k=1;k<=j;k++) sA(j,k) -= conj(se(j))*sA(i,k)+conj(sA(i,j))*se(k);
    }
    L270:
    for (k=1;k<=l;k++) sA(i,k) *= scale;

    se(l) = std::complex<T>(real(se(l)),-imag(stau(l)));
    L290:
    hh = sd(i);
    sd(i) = real(sA(i,i));
    sA(i,i) = std::complex<T>(hh,scale * sqrt(h));
  };//L300:
	evals=real(ee);
	for (int i=1;i<=n;i++) stau(i)=std::complex<T>(real(stau(i)),imag(se(i)));
  return;
}


typedef TYPE Type;
typedef MTYPE<std::complex<Type> > MType;
template void htridi<Type,MType>(MType& A, Vector<Type>& d ,Vector<Type>& evals,
	    Vector<Type>& e2,Vector<std::complex<Type> >& tau);
