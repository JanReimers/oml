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
//  SMatrix             MatrixShape   Sym     Real.
//  HMatrix             MatrixShape   Sym     Real.
//  GetColumn(Matrix)   VectorShape   Full    Abstract.
//  Transpose(Matrix)   MatrixShape   Full    Abstract.
//  OuterProduct(V1,V2) MatrixShape   Sym     Abstract.
//

enum Shape {ArrayShape,VectorShape,MatrixShape};
enum Store {Full,Symmetric,Diagonal,Sparse};
enum Data  {Real,Abstract};

// Define how different storage layouts and Data type interact
template <Store S1, Store S2> class ReturnStore{};
template <Store S> class ReturnStore<S   ,S   >{public: static const Store RetType=S;};
template <Store S> class ReturnStore<Full,S   >{public: static const Store RetType=Full;};
template <Store S> class ReturnStore<S   ,Full>{public: static const Store RetType=Full;};
template <       > class ReturnStore<Full,Full>{public: static const Store RetType=Full;};
template <>        class ReturnStore<Symmetric,Diagonal >{public: static const Store RetType=Symmetric;};
template <>        class ReturnStore<Diagonal ,Symmetric>{public: static const Store RetType=Symmetric;};

template <Data D1, Data D2> class ReturnData;
template <Data D> class ReturnData<D       ,D       > {public: static const Data RetType=D;};
template <Data D> class ReturnData<D       ,Abstract> {public: static const Data RetType=Abstract;};
template <Data D> class ReturnData<Abstract,D       > {public: static const Data RetType=Abstract;};
template <>       class ReturnData<Abstract,Abstract> {public: static const Data RetType=Abstract;};

#endif // _shape_h_
