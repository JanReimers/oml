module;
#include <vector>

export module omlNumericalRecipes;
// import omlNumericalRecipes:Imp;
import oml.Vector;
import oml.SMatrix;
import oml.Matrix;

export{
//----------------------------------------------------------------------------
/*! \file numeric.h
  \brief Numerical methods as template functions.
*/
//! Cholsky decomposition of a symmetric matrix, that is \b not stored symmetrically.
template <class T> void Cholsky          (Matrix <T>&                         );
//! Cholsky decomposition of a symmetric matrix, that is stored symmetrically.
template <class T> void Cholsky          (SMatrix<T>&                         );
//! Sort columns of eigenstd::vectors based on the order of eigen values.
template <class T,class M> void EigenSort        (M&, Vector<T>& EigenValues);
//! Solve Ax=B using GJ elimintaion. x is returned in A.
template <class T> void GaussJordanSolver(Matrix <T>& A, Matrix<T>& B          );
//! Invert an upper triangular matrix.
template <class T> void InvertTriangular (Matrix <T>&                         );
//! Invert an upper triangular matrix stored symmetrically.
template <class T> void InvertTriangular (SMatrix<T>&                         );
//! Solve tri-diagonal system TDx=A. x is returned in A.
template <class T, class M> void TriDiagonal      (M& A, Vector<T>& Diagonal, Vector<T>& OffDiagonal);
//! Diagonalize matrix A. Eigen std::vectors returned in A, eigen value returned as a new Vector.
template <class T> Vector<T> Diagonalize(Matrix<T>& A);

//! LU back substitution
template <class T> void LUBackSub(const Matrix<T>& A, Vector<T>& B,const std::vector<index_t>& Index);
//! LU back substitution
template <class T> void LUBackSub(const Matrix<T>& A, Matrix<T>& B,const std::vector<index_t>& Index);
//! LU decomposition of A=L*U.
template <class T> bool LUDecomp(Matrix<T>& A, std::vector<index_t>& Index ,T& d);
//! QL decomposition.  Used for eigen problems.
template <class T,class M> void QLDecomp( M& A, Vector<T>& Diagonal, Vector<T>& OffDiagonal);
//! Single value back substitution.
template <class T> void SVBackSub(const Matrix<T>& U, const Vector<T>& W,const Matrix<T>& V, const Vector<T>& B, Vector<T>& X);
//! Single value decomposition.
template <class T, class M> void SVDecomp (M& A,Vector<T>& W,M& V);

//! Invert a square matrix (not neccesarily symmetric).
template <class T> Matrix <T> InvertSquare   (const Matrix <T>&);
//! Invert a symmtric matrix with full sotrage.
template <class T> Matrix <T> InvertSymmetric(const Matrix <T>&);
//! Invert a symmtric matrix with symmetric sotrage.
template <class T> SMatrix<T> InvertSymmetric(const SMatrix<T>&);

//! Solve symmtric tri diagonal system Tu=r where T is a symmetric tridiagonal matrix. Last element of OffDiagonal is not used.
//Array<double> SolveTriDiag(const Array<double> & OffDiagonal,
//			   const Array<double> & Diagonal,
//			   const Array<double> & r);

//! Solve symmtric tri diagonal system Tu=r where T is a tridiagonal matrix. Last element of OffDiagonal is not used.
Vector<double> SolveTriDiag(const Vector<double> & OffDiagonal,
			    const Vector<double> & Diagonal,
			    const Vector<double> & r);

template <class T> Vector<T> Diagonalize(Matrix<T>& m)
{
  assert(m.GetRowLimits()==m.GetColLimits());
  assert(m.GetRowLimits()==m.GetColLimits());
  assert(IsSymmetric(m));

  Vector<T> EigenValues(m.GetRowLimits());
  Vector<T> OffDiagonal(m.GetRowLimits());

  TriDiagonal(m,EigenValues,OffDiagonal);
  QLDecomp   (m,EigenValues,OffDiagonal);
  EigenSort  (m,EigenValues);
  return EigenValues;
}

template <class T> std::tuple<Matrix<T>,Vector<T> > Diagonalize(const SMatrix<T>& s)
{
  Vector<T> EigenValues(s.GetRowLimits());
  Vector<T> OffDiagonal(s.GetRowLimits());

  Matrix<T>  m(s);
  TriDiagonal(m,EigenValues,OffDiagonal);
  QLDecomp   (m,EigenValues,OffDiagonal);
  EigenSort  (m,EigenValues);
  return std::make_tuple(m,EigenValues);
}

template <class T> std::tuple<Matrix<T>,Vector<T>,Matrix<T> > SVD(const SMatrix<T>& A)
{
    Matrix<T> U(A),V(A.GetLimits());
    Vector<T> s(A.GetRowLimits());
    SVDecomp (U,s,V);
    return std::make_tuple(U,s,V);
}



} // export block


// using Type=double;
// template void Cholsky<Type>(Matrix<Type>& A);
// template void EigenSort(Matrix<Type>&, Vector<Type>&);
// template  Matrix<Type> InvertSquare   (const  Matrix<Type>&);
// template  Matrix<Type> InvertSymmetric(const  Matrix<Type>&);
// template SMatrix<Type> InvertSymmetric(const SMatrix<Type>&);
// template void InvertTriangular(Matrix<Type>&);
// template void InvertTriangular(SMatrix<Type>&);
// template void LUBackSub(const Matrix<Type>&, Vector<Type>&,const std::vector<index_t>&);
// template void LUBackSub(const Matrix<Type>&, Matrix<Type>&, const std::vector<index_t>&);
// template bool LUDecomp(Matrix<Type>&,std::vector<index_t>&,Type&);
// template void QLDecomp(Matrix<Type>&,Vector<Type>&,Vector<Type>&);
// template void SVBackSub(const Matrix<Type>&,const Vector<Type>&,const Matrix<Type>&,const Vector<Type>&, Vector<Type>&);
// template void SVDecomp<Type,Matrix<Type>>(Matrix<Type>&,Vector<Type>&,Matrix<Type>&);
// template void TriDiagonal(Matrix<Type>&,Vector<Type>&,Vector<Type>&);