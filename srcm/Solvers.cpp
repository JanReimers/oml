module;
#include <tuple>
#include <complex>

export module oml.Solvers;

import oml.Matrix;
import oml.DiagonalMatrix;
import oml.Vector;

//gcc15.0 bug it won't instantiate this by itsef in release mode (-O2).
template class std::tuple<Matrix<std::complex<double>>,Vector<std::complex<double>>>;
export namespace oml
{
//
//  Abstract interface for eigen solvers
//
template <class T> class EigenSolver
{
protected:
    typedef std::complex<double> dcmplx;
    typedef Vector <double> VectorR;
    typedef Vector <dcmplx> VectorC;
    typedef Matrix<T>       MatrixT;
    typedef Matrix<dcmplx>  MatrixC;
    typedef std::tuple<MatrixT,VectorR> UdType;
    typedef std::tuple<MatrixC,VectorC> UdTypeN;
public:
    virtual ~EigenSolver() {};
    virtual void Reset()=0;  //Clear out guesses from previous iterations

    virtual UdType  Solve              (const MatrixT&,double eps, int NumEigenValues)=0;
    virtual UdType  SolveAll           (const MatrixT&,double eps                    )=0;
    virtual UdTypeN SolveRightNonSym   (const MatrixT&,double eps, int NumEigenValues)=0;
    virtual UdTypeN SolveLeft_NonSym   (const MatrixT& A,double eps, int NumEigenValues)
    {
        MatrixT Adagger=Transpose(conj(A));
        auto [U,d]=SolveRightNonSym(Adagger,eps,NumEigenValues);
        return std::make_tuple(conj(U),conj(d));
    }

    virtual UdTypeN SolveAllRightNonSym(const MatrixT&,double eps                    )=0;

};

//
//  Abstract interface for Linear solvers  Ax=b.
//  The combinatorics of Matrix shapes and data types can be overwhelming
//  So we just add things as we need them.
//
template <class T> class LinearSolver
{
protected:
    typedef Matrix<T> MatrixT;
    typedef Vector<T> VectorT;
public:
    virtual ~LinearSolver() {};
//
//  Non packed triangular systems
//      solve A*x=b
//
    virtual VectorT SolveUpperTri(const MatrixT& A,const VectorT& b)=0; // A is upper triangular
    virtual VectorT SolveLowerTri(const MatrixT& A,const VectorT& b)=0; // A is lower triangular
//
//      solve b=x*A
//
//    virtual VectorT SolveUpperTri(const VectorT& b,const MatrixT& A); // A is upper triangular
//    virtual VectorT SolveLowerTri(const VectorT& b,const MatrixT& A); // A is lower triangular
//
//  Full solver
//
    virtual VectorT Solve(const MatrixT& A,const VectorT& b)=0; // A*x=b
//    virtual VectorT Solve(const VectorT& b,const MatrixT& A);   // b=x*A
};

//
//  Abstract interface for SVD solvers
//
template <class T> class SVDSolver
{
protected:
    typedef         Matrix<T>               MatrixT;
    typedef DiagonalMatrix<double> DiagonalMatrixRT;
    typedef std::tuple<MatrixT,DiagonalMatrixRT,MatrixT> UsVType;
public:
    virtual ~SVDSolver() {};

    virtual UsVType  Solve   (const MatrixT&, size_t NumSingluarValues)=0;
    virtual UsVType  SolveAll(const MatrixT&                       )=0;

};



} // namespace
