//---------------------------------------------------------------------------------
//
//  Make template instances
//
#include "src/matrix.cpp"
#include "src/vector.cpp"
#define Type double

template class DMatrix<Type>;
template class  Vector<Type>;

template void SetLimits     (Indexable<Type,DMatrix<Type>,Full,Real,MatrixShape>&,const MatLimits&,bool);
template void ReIndexRows   (Indexable<Type,DMatrix<Type>,Full,Real,MatrixShape>&,const std::vector<index_t>& index);
template void ReIndexColumns(Indexable<Type,DMatrix<Type>,Full,Real,MatrixShape>&,const std::vector<index_t>& index);
template void SwapRows      (Indexable<Type,DMatrix<Type>,Full,Real,MatrixShape>& m,index_t i,index_t j);
template void SwapColumns   (Indexable<Type,DMatrix<Type>,Full,Real,MatrixShape>& m,index_t i,index_t j);
template void SubMatrix     (Indexable<Type,DMatrix<Type>,Full,Real,MatrixShape>& m,const Indexable<Type,DMatrix<Type>,Full,Real,MatrixShape>& old);





