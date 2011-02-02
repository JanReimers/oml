#include "oml/vector.h"
#include "oml/matrix.h"
#include "oml/cnumeric.h"
#include <cassert>
#include <iostream>

// This code comes from netlib.
// Modifications for C++ and oml containers Copyright (1994-2003), Jan N. Reimers

//---------------------------------------------------------------------
//
//     this subroutine calls the recommended sequence of
//     subroutines from the eigensystem subroutine package (eispack)
//     to find the eigenvalues and eigenvectors (if desired)
//     of a complex hermitian matrix.
//
void ch(Matrix<std::complex<double> >& A,Vector<double>& w ,bool matz,int ierr)
{
  assert(A.GetNumRows()==A.GetNumRows());
  int n=A.GetNumRows();
  
  Vector<double> fv1(n),fv2(n);
  Vector<std::complex<double> > fm1(n);
  
  htridi(A,w,fv1,fv2,fm1);
   if (!matz)
     tqlrat(w,fv2,ierr); //find eigenvalues only
   else
   {	
     Matrix<std::complex<double> > Z(A.GetLimits());
     Unit(Z);
     tql2(w,fv1,Z,ierr);
     htribk(A,fm1,Z);
     A=Z;
   }
  
   return;
}
