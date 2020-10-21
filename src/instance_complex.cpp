//---------------------------------------------------------------------------------
//
//  Make template instance
//
#include "src/matrix.cpp"
#include "src/vector.cpp"
#define Type std::complex<double>

template class Matrix<Type>;
template class  Vector<Type>;

template void SetLimits     (Indexable<Type,Matrix<Type>,Full,Real,MatrixShape>&,const MatLimits&,bool);
template void ReIndexRows   (Indexable<Type,Matrix<Type>,Full,Real,MatrixShape>&,const std::vector<index_t>& index);
template void ReIndexColumns(Indexable<Type,Matrix<Type>,Full,Real,MatrixShape>&,const std::vector<index_t>& index);
template void SwapRows      (Indexable<Type,Matrix<Type>,Full,Real,MatrixShape>& m,index_t i,index_t j);
template void SwapColumns   (Indexable<Type,Matrix<Type>,Full,Real,MatrixShape>& m,index_t i,index_t j);
template void SubMatrix     (Indexable<Type,Matrix<Type>,Full,Real,MatrixShape>& m,const Indexable<Type,Matrix<Type>,Full,Real,MatrixShape>& old);

