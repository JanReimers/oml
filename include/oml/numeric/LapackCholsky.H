#ifndef _LapackCholsky_H_
#define _LapackCholsky_H_

template <class T> class SMatrix;
template <class T> class Matrix;

namespace oml {

// This should really return an upper triangular matrix, but we don't support that storage type yet.
// The input matrices must be positive definite.

 Matrix<double> LapackCholsky(const SMatrix<double>&);
SMatrix<double> LapackInvertSymmetric(const SMatrix<double>&); 


};

#endif // _LapackCholsky_H_
