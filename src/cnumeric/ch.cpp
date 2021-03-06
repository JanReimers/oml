#include "oml/vector.h"
#include "oml/matrix.h"
#include "oml/cnumeric.h"
#include <cassert>
#include <iostream>
#include <complex>

#ifndef TYPE
#error "ch.cc TYPE not defined"
#endif

#ifndef MTYPE
#error "ch.cc MTYPE not defined"
#endif


// This code comes from netlib.
// Modifications for C++ and oml containers Copyright (1994-2003), Jan N. Reimers

//---------------------------------------------------------------------
//
//     this subroutine calls the recommended sequence of
//     subroutines from the eigensystem subroutine package (eispack)
//     to find the eigenvalues and eigenvectors (if desired)
//     of a complex hermitian matrix.
//
template <class T, class M>
void ch(M& A,Vector<T>& w ,bool matz,int& ierr)
{
  assert(A.GetNumRows()==A.GetNumRows());
  index_t n=A.GetNumRows();

  Vector<T> fv1(n),fv2(n);
  Vector<std::complex<T> > fm1(n);

  htridi(A,w,fv1,fv2,fm1);
   if (!matz)
     tqlrat(w,fv2,ierr); //find eigenvalues only
   else
   {
     M Z(A.GetLimits());
     Unit(Z);
     tql2(w,fv1,Z,ierr);
     htribk(A,fm1,Z);
     A=Z;
   }

   return;
}

typedef TYPE Type;
typedef MTYPE<std::complex<Type> > MType;
template void ch<Type,MType>(MType&,Vector<Type>&,bool,int&);


