#ifndef _LapackCholsky_H_
#define _LapackCholsky_H_

template <class T> class SMatrix;
template <class T> class Matrix;

namespace oml {

// This should really return an upper triangular matrix, but we dont
Matrix<double> LapackCholsky(const SMatrix<double>&);

};

#endif // _LapackCholsky_H_
