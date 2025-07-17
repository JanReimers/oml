// File: UT-Eigen-double.cc  Time the eigen std::vector routine for Matrix<double>

// Copyright (1994-2003), Jan N. Reimers

#include "oml/UnitTests/UnitTest.H"
#include "oml/stopw.h"
#include "oml/vector.h"
#include "oml/matrix.h"
#include "oml/numeric.h"
#include <iostream>
#include <cmath>

int TestEigenDouble()
{
  const char* Class="Eigen tester";
  StartClass(Class);

  int N=500;

  Matrix<double> A(N,N);
  Vector<double> E(N);
  StopWatch sw;

  for (int i=A.GetRowLow();i<=A.GetRowHigh();i++)
    for (int j=i;j<=A.GetColHigh();j++)
    {
      double del=(double)(i-j)/N;
      A(j,i)=A(i,j)=exp(-del*del);
    }
//  A.Random(1.0);
//  A.Symmetric();

  if(!IsSymmetric(A))
  {
    std::cerr << "A is nont symmetric!!" << std::endl;
    return -1;
  }

  sw.Start();
  E=Diagonalize(A);
  sw.Stop();

  std::cout << N << "X" << N << " matrix diagonalized in "<< sw.GetTime() <<" sec" << std::endl;
  return 0;
}
