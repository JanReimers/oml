#ifndef QRSOLVER_H
#define QRSOLVER_H

#include <tuple>

template <class T> class Matrix;
//
//  Abstract interface for QR solvers
//
template <class T> class QRSolver
{
protected:
    typedef Matrix<T>                   MatrixT;
    typedef std::tuple<MatrixT,MatrixT> QRType;
public:
    virtual ~QRSolver() {};

    virtual QRType  SolveThinQR(const MatrixT&)=0; // returns [Q,R]
    virtual QRType  SolveThinRQ(const MatrixT&)=0; // mirror image version returns [R,Q]
    virtual QRType  SolveThinQL(const MatrixT&)=0; // returns [Q,L]
    virtual QRType  SolveThinLQ(const MatrixT&)=0; // mirror image version returns [L,Q]

    virtual QRType  SolveRankRevealingQR(const MatrixT&, double eps)=0; // returns [Q,R]
    virtual QRType  SolveRankRevealingRQ(const MatrixT&, double eps)=0; // returns [Q,L]
    virtual QRType  SolveRankRevealingQL(const MatrixT&, double eps)=0; // returns [R,Q]
    virtual QRType  SolveRankRevealingLQ(const MatrixT&, double eps)=0; // returns [L,Q]

    virtual QRType  SolveThinQR(const MatrixT& A,double eps)
    {
        if (eps>=0.0)
            return SolveRankRevealingQR(A,eps);
        else
            return SolveThinQR(A);
    }
    virtual QRType  SolveThinRQ(const MatrixT& A,double eps)
    {
        if (eps>=0.0)
            return SolveRankRevealingRQ(A,eps);
        else
            return SolveThinRQ(A);
    }
    virtual QRType  SolveThinQL(const MatrixT& A,double eps)
    {
        if (eps>=0.0)
            return SolveRankRevealingQL(A,eps);
        else
            return SolveThinQL(A);
    }
    virtual QRType  SolveThinLQ(const MatrixT& A,double eps)
    {
        if (eps>=0.0)
            return SolveRankRevealingLQ(A,eps);
        else
            return SolveThinLQ(A);
    }
};

#endif // QRSOLVER_H
