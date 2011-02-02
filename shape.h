// File: shape.h  Abstract interface for communicating container shapes.
#ifndef _shape_h_
#define _shape_h_

// Copyright (1994-2003), Jan N. Reimers

//--------------------------------------------------
//
//  Containers will have the following traits:
//                       Shape        Store   Data
//  ------------------------------------------------
//  Array               ArrayShape    Full    Real.
//  Vector              VectorShape   Full    Real.
//  Matrix              MatrixShape   Full    Real.
//  DMatrix             MatrixShape   Full    Real.
//  SMatrix             MatrixShape   Sym     Real.
//  HMatrix             MatrixShape   Sym     Real.
//  GetColumn(Matrix)   VectorShape   Full    Abstract.
//  Transpose(Matrix)   MatrixShape   Full    Abstract.
//  OuterProduct(V1,V2) MatrixShape   Sym     Abstract.
//

enum Shape {ArrayShape,VectorShape,MatrixShape};
enum Store {Full,Symmetric};
enum Data  {Real,Abstract};

#endif // _shape_h_
