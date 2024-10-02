//#include "NumericalMethods/SparseSVDSolver.H"
//#include "NumericalMethods/SparseEigenSolver.H"
#include "oml/numeric/EigenSolver.H"
//#include "NumericalMethods/LinearSolver.H"
//#include "Containers/SparseMatrix.H"
#include "oml/matrix.h"
#include "oml/vector.h"
#include "oml/fakedouble.h"

namespace oml {

template <class T> typename EigenSolver<T>::UdTypeN
EigenSolver<T>::SolveLeft_NonSym(const MatrixT& A,double eps, int NumEigenValues)
{
    MatrixT Adagger=Transpose(conj(A));
    auto [U,d]=SolveRightNonSym(Adagger,eps,NumEigenValues);
    return std::make_tuple(conj(U),conj(d));
}

typedef std::complex<double> dcmplx;

template class EigenSolver<double>;
template class EigenSolver<dcmplx>;

} //namepsace oml
