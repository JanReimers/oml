#include "Epsilons.H"
#include <iostream>
#include <sstream>
#include <iomanip>

namespace oml
{

Epsilons::Epsilons(double default_eps)
    : itsDelatEnergy1Epsilon        (default_eps)
    , itsDelatEnergy2Epsilon        (default_eps)
    , itsDelatNormEpsilon           (default_eps)
    , itsEigenSolverEpsilon         (default_eps)
    , itsNormalizationEpsilon       (default_eps)
    , itsMPSCompressEpsilon         (default_eps)
    , itsMPOCompressEpsilon         (default_eps)
    , itsSparseMatrixEpsilon        (default_eps)
    , itsDeltaLambdaEpsilon         (default_eps)
{}


std::ostream& operator<<(std::ostream& os,const Epsilons& eps)
{
    os.precision(1);
    const static int width=3;
    os
    << std::setw(width) << std::scientific << eps.itsDelatEnergy1Epsilon  << " "
    << std::setw(width) << std::scientific << eps.itsDelatEnergy2Epsilon  << " "
    << std::setw(width) << std::scientific << eps.itsDelatNormEpsilon     << " "
    << std::setw(width) << std::scientific << eps.itsEigenSolverEpsilon   << " "
    << std::setw(width) << std::scientific << eps.itsNormalizationEpsilon << " "
    << std::setw(width) << std::scientific << eps.itsMPSCompressEpsilon   << " "
    << std::setw(width) << std::scientific << eps.itsMPOCompressEpsilon   << " "
    << std::setw(width) << std::scientific << eps.itsSparseMatrixEpsilon  << " "
    << std::setw(width) << std::scientific << eps.itsDeltaLambdaEpsilon   << " "
    ;
    return os;
}

std::string Epsilons::Header()
{
    std::ostringstream os;
    os << " <E>    <E^2>-<E>^2 1-<|>  Eigen  Norm    SVD1  SVD2   Sparse dLambda" << std::ends;
    return os.str();
}

}
