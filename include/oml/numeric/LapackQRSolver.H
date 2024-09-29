#ifndef LAPACKQRSOLVER_H
#define LAPACKQRSOLVER_H

#include "NumericalMethods/QRSolver.H"

//
//  Abstract interface for QR solvers
//
template <class T> class LapackQRSolver
    : public virtual QRSolver<T>
{
protected:
    using MatrixT=typename QRSolver<T>::MatrixT;
    using QRType =typename QRSolver<T>::QRType;
public:
    virtual ~LapackQRSolver() {};

//    virtual QRType  Solve    (const MatrixT&);
    virtual QRType  SolveThinQR(const MatrixT&); // returns [Q,R]
    virtual QRType  SolveThinRQ(const MatrixT&); // mirror image version returns [R,Q]
    virtual QRType  SolveThinQL(const MatrixT&); // returns [Q,L]
    virtual QRType  SolveThinLQ(const MatrixT&); // mirror image version returns [L,Q]
//
//  Rank revealing QR removes zero rows (within eps) from the R matrix and the corresponding columns
//  of Q.
//
    virtual QRType  SolveRankRevealingQR(const MatrixT&, double eps); // returns [Q,R]
    virtual QRType  SolveRankRevealingRQ(const MatrixT&, double eps); // returns [Q,L]
    virtual QRType  SolveRankRevealingQL(const MatrixT&, double eps); // returns [Q,R]
    virtual QRType  SolveRankRevealingLQ(const MatrixT&, double eps); // returns [L,Q]

};

#endif // LAPACKQRSOLVER_H
