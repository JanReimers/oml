//#include "NumericalMethods/SparseSVDSolver.H"
//#include "NumericalMethods/SparseEigenSolver.H"
#include "oml/numeric/EigenSolver.H"
//#include "NumericalMethods/LinearSolver.H"
//#include "Containers/SparseMatrix.H"
#include "oml/matrix.h"
#include "oml/vector.h"
#include "oml/fakedouble.h"

template <class T> typename EigenSolver<T>::UdTypeN
EigenSolver<T>::SolveLeft_NonSym(const MatrixT& A,double eps, int NumEigenValues)
{
    MatrixT Adagger=Transpose(conj(A));
    auto [U,d]=SolveRightNonSym(Adagger,eps,NumEigenValues);
    return std::make_tuple(conj(U),conj(d));
}

//template <class T> typename SparseEigenSolver<T>::UdTypeN
//SparseEigenSolver<T>::SolveLeft_NonSym(const SparseMatrixT& A,double eps, int NumEigenValues)
//{
//    SparseMatrixT Adagger=~A;
//    auto [U,d]=SolveRightNonSym(Adagger,eps,NumEigenValues);
//    return std::make_tuple(conj(U),conj(d));
//}
//
//template <class T> typename LinearSolver<T>::VectorT
//LinearSolver<T>::SolveUpperTri(const VectorT& b,const MatrixT& A)
//{
//    return SolveLowerTri(Transpose(A),b);
//}
//template <class T> typename LinearSolver<T>::VectorT
//LinearSolver<T>::SolveLowerTri(const VectorT& b,const MatrixT& A)
//{
//    return SolveUpperTri(Transpose(A),b);
//}
//template <class T> typename LinearSolver<T>::VectorT
//LinearSolver<T>::Solve(const VectorT& b,const MatrixT& A)
//{
//    return Solve(Transpose(A),b);
//}
//
//  static variable. Kludge for getting the matrix into the Mat*Vec routines.
//
//template<class T> const SparseMatrix<T>* SparseEigenSolver<T>::theSparseMatrix = 0;
//template<class T> const       Matrix<T>* SparseEigenSolver<T>::theDenseMatrix = 0;
//template<class T> const SparseMatrix<T>* SparseSVDSolver<T>::theSparseMatrix = 0;
//template<class T> const       Matrix<T>* SparseSVDSolver<T>::theDenseMatrix = 0;
//
//  Used to pass client pointers in mat*vec functions
//
//template <class T> const SparseEigenSolverClient<T>* SparseEigenSolverClient<T>::theClient;
//template <class T> const SparseSVDSolverClient<T>* SparseSVDSolverClient<T>::theClient;

typedef std::complex<double> dcmplx;

//template class SparseSVDSolver<double>;
//template class SparseSVDSolver<dcmplx>;
template class EigenSolver<double>;
template class EigenSolver<dcmplx>;
//template class SparseEigenSolver<double>;
//template class SparseEigenSolver<dcmplx>;
//template class SparseSVDSolverClient<double>;
//template class SparseSVDSolverClient<dcmplx>;
//template class SparseEigenSolverClient<double>;
//template class SparseEigenSolverClient<dcmplx>;
//template class LinearSolver<double>;
//template class LinearSolver<dcmplx>;

/*
//#include "TensorNetworks/Enums.H"
//
//  Search right to left for the last non zero element, and return it's index.
//
index_t FindLast(const Vector<double>& v,double eps)
{
    VecLimits l=v.GetLimits();
    index_t last_index=l.Low;
    for (index_t i=l.High;i>l.Low;i--)
         if (fabs(v(i))>eps)
         {
             last_index=i;
             break;
         }
     return last_index;
}
//
//  Search left to right for the first non zero element, and return it's index.
//
index_t FindFirst(const Vector<double>& v,double eps)
{
    VecLimits l=v.GetLimits();
    index_t first_index=l.High;
    for (index_t i=l.Low;i<l.High;i++)
         if (fabs(v(i))>eps)
         {
             first_index=i;
             break;
         }
     return first_index;
}

//
//  Create a re-indexing array to bring the matrix into Upper/Lower triangular form
//
std::vector<index_t>  FindRowReIndex(TensorNetworks::TriType ul,const Matrix<double>& UL, double eps)
{
    MatLimits l=UL.GetLimits();
    //
    //  Search each row for the first/last non-zero element.
    //        col     row
    std::multimap<index_t,index_t> index;
    for (index_t ir:UL.rows())
    {
        if      (ul==Lower) index.insert({FindLast (UL.GetRow(ir),eps),ir});
        else if (ul==Upper) index.insert({FindFirst(UL.GetRow(ir),eps),ir});
        else assert(false);
    }
    //
    //  Build up the re-index array from <col,row> map.
    //
    std::vector<index_t> reindex;
    for (const auto& [key, value] : index)
        reindex.push_back(value-l.Row.Low); //Row numbers need to be zero based for the oml ReIndex routines.

    return reindex;
}
std::vector<index_t>  FindColReIndex(TensorNetworks::TriType ul,const Matrix<double>& U, double eps)
{
    MatLimits l=U.GetLimits();
    //
    //  Search each col for the first/last non-zero element.
    //        row     col
    std::multimap<index_t,index_t> index;
    for (index_t ic:U.cols())
    {
        if      (ul==Lower) index.insert({FindFirst(U.GetColumn(ic),eps),ic});
        else if (ul==Upper) index.insert({FindLast (U.GetColumn(ic),eps),ic});
        else assert(false);
    }
    //
    //  Build up the re-index array from <row,col> map.
    //
    std::vector<index_t> reindex;
    for (const auto& [key, value] : index)
        reindex.push_back(value-l.Col.Low); //Col numbers need to be zero based for the oml ReIndex routines.

    return reindex;
}

*/

