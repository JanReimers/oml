#ifndef EPSILONS_H
#define EPSILONS_H

#include <iosfwd>

namespace TensorNetworks
{
//
//  Simple struct for passing many epsilons down through the call hierarchy
//
struct Epsilons
{
    Epsilons(double default_eps=1e-12);
    double itsDelatEnergy1Epsilon;      // assume converged when <E_n> - <E_n-1> < eps
    double itsDelatEnergy2Epsilon;      // assume converged when <E^2> - <E>^2 < eps this tests convergence in D
    double itsDelatNormEpsilon;         // For Psi optimization iterations after SVD compression.
    double itsEigenSolverEpsilon;       // Decides convergence of whatever eigen solver is used.
    double itsNormalizationEpsilon;     // Is normalized if max|A-I|<eps
    double itsMPSCompressEpsilon;       // s is deemed to be zero if s<eps
    double itsMPOCompressEpsilon;       // s is deemed to be zero if s<eps
    double itsSparseMatrixEpsilon;      // matrix element A_ij is deemed to be zero if A_ij <eps
    double itsDeltaLambdaEpsilon;       // assume convergence when delta lambdaA or B < eps

    friend std::ostream& operator<<(std::ostream&,const Epsilons&);
    static std::string Header();
};

}
#endif // EPSILONS_H
