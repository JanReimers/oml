#include "Tests.H"

#include "Epsilons.H"
#include "LinearAlgebraTests.H"
//#include "NumericalMethods/PrimeSVDSolver.H"
//#include "NumericalMethods/PrimeEigenSolver.H"
#include "oml/numeric/LapackSVDSolver.H"
#include "oml/numeric/LapackEigenSolver.H"
//#include "NumericalMethods/ArpackEigenSolver.H"
//#include "NumericalMethods/LapackQRSolver.H"
//#include "NumericalMethods/LapackLinearSolver.H"
//#include "Containers/SparseMatrix.H"
#include "oml/numeric.h"
#include "oml/cnumeric.h"
#include "oml/diagonalmatrix.h"

class LinearAlgebraTests : public ::testing::Test
{
public:
    typedef Matrix<dcmplx> MatrixCT;
    typedef Matrix<double> MatrixRT;
    typedef Vector <double> VectorRT;
    typedef DiagonalMatrix<double> DiagonalMatrixRT;
//    typedef SparseMatrix<dcmplx> SparseMatrixCT;
//    typedef SparseMatrix<double> SparseMatrixRT;

    LinearAlgebraTests()
    : itsEps()
    , eps(1.0e-13)
#ifdef DEBUG
    , Nsvd   (30)
    , Neigen (30)
    , Nqr    (30)
    , Nlinear(30)
    , Ncholsky(30)
#else
    , Nsvd   (100)
    , Neigen (100)
    , Nqr    (200)
    , Nlinear(200)
    , Ncholsky(100)
#endif
    , svdDensity(0.2)
    , eigenDensity(0.1)
    {
        StreamableObject::SetToPretty();
    }
    void SetupC(int M,int N)
    {
        int mn=Min(M,N);
        itsAC.SetLimits(M,N);
        itsIC.SetLimits(mn,mn);
        FillRandom(itsAC);
        Unit(itsIC);
    }
    void SetupH(int N)
    {
        MatrixCT A(N,N);
        itsWR.SetLimits(N);
        itsIC.SetLimits(N,N);
        FillRandom(A);
        itsAC=A+Transpose(conj(A)); //Make it hermitian
        Unit(itsIC);
    }

//    SparseMatrixRT itsARs;
//    SparseMatrixCT itsACs;
    MatrixCT  itsAC,itsIC;
    MatrixRT  itsAR,itsIR;
    VectorRT  itsWR; //Eigen values

    TensorNetworks::Epsilons  itsEps;
    double eps;
    int Nsvd,Neigen,Nqr,Nlinear,Ncholsky;
    double svdDensity;
    double eigenDensity;

};

//
//TEST_F(LinearAlgebraTests,SparseMatrix)
//{
//    Matrix<double> d(5,6);
//    FillRandom(d);
//    SparseMatrix<double> m=d;
////    cout << m << endl;
//}
//
//TEST_F(LinearAlgebraTests,Primme_SVDSolverDenseReal)
//{
//    SVDTester<double>(new PrimeSVDSolver<double>(),Nsvd  ,Nsvd  ).RunTests();
//    SVDTester<double>(new PrimeSVDSolver<double>(),Nsvd/2,Nsvd  ).RunTests();
//    SVDTester<double>(new PrimeSVDSolver<double>(),Nsvd  ,Nsvd/2).RunTests();
//}
//
//TEST_F(LinearAlgebraTests,Primme_SVDSolverDenseComplex)
//{
//    SVDTester<dcmplx>(new PrimeSVDSolver<dcmplx>(),Nsvd  ,Nsvd  ).RunTests();
//    SVDTester<dcmplx>(new PrimeSVDSolver<dcmplx>(),Nsvd/2,Nsvd  ).RunTests();
//    SVDTester<dcmplx>(new PrimeSVDSolver<dcmplx>(),Nsvd  ,Nsvd/2).RunTests();
//}
//
//TEST_F(LinearAlgebraTests,Primme_SVDSolverSparseReal)
//{
//    SparseSVDTester<double>(new PrimeSVDSolver<double>(),Nsvd  ,Nsvd  ,svdDensity).RunTests();
//    SparseSVDTester<double>(new PrimeSVDSolver<double>(),Nsvd/2,Nsvd  ,svdDensity).RunTests();
//    SparseSVDTester<double>(new PrimeSVDSolver<double>(),Nsvd  ,Nsvd/2,svdDensity).RunTests();
//}
//
//TEST_F(LinearAlgebraTests,Primme_SVDSolverSparseComplex)
//{
//    SparseSVDTester<dcmplx>(new PrimeSVDSolver<dcmplx>(),Nsvd  ,Nsvd  ,svdDensity).RunTests();
//    SparseSVDTester<dcmplx>(new PrimeSVDSolver<dcmplx>(),Nsvd/2,Nsvd  ,svdDensity).RunTests();
//    SparseSVDTester<dcmplx>(new PrimeSVDSolver<dcmplx>(),Nsvd  ,Nsvd/2,svdDensity).RunTests();
//}
//
TEST_F(LinearAlgebraTests,Lapack_SVDSolverDenseReal)
{
    SVDTester<double>(new oml::LapackSVDSolver<double>(),Nsvd  ,Nsvd  ).RunTests();
    SVDTester<double>(new oml::LapackSVDSolver<double>(),Nsvd/2,Nsvd  ).RunTests();
    SVDTester<double>(new oml::LapackSVDSolver<double>(),Nsvd  ,Nsvd/2).RunTests();
}
TEST_F(LinearAlgebraTests,Lapack_SVDSolverDenseComplex)
{
    SVDTester<dcmplx>(new oml::LapackSVDSolver<dcmplx>(),Nsvd  ,Nsvd  ).RunTests();
    SVDTester<dcmplx>(new oml::LapackSVDSolver<dcmplx>(),Nsvd/2,Nsvd  ).RunTests();
    SVDTester<dcmplx>(new oml::LapackSVDSolver<dcmplx>(),Nsvd  ,Nsvd/2).RunTests();
}


//TEST_F(LinearAlgebraTests,Primme_EigenSolverDenseReal)
//{
//    EigenTester<double>(new PrimeEigenSolver<double>(),Neigen).RunTests(true);
//}
//
//TEST_F(LinearAlgebraTests,Primme_EigenSolverDenseComplex)
//{
//    EigenTester<dcmplx>(new PrimeEigenSolver<dcmplx>(),Neigen).RunTests(true);
//}

//TEST_F(LinearAlgebraTests,Primme_EigenSolverSparseReal)
//{
//    SparseEigenTester<double>(new PrimeEigenSolver<double>(),Neigen,eigenDensity).RunTests(true);
//}
//
//TEST_F(LinearAlgebraTests,Primme_EigenSolverSparseComplex)
//{
//    SparseEigenTester<dcmplx>(new PrimeEigenSolver<dcmplx>(),Neigen,eigenDensity).RunTests(true);
//}
//

TEST_F(LinearAlgebraTests,Lapack_EigenSolverDenseReal)
{
    EigenTester<double>(new oml::LapackEigenSolver<double>(),Neigen).RunTests();
}
TEST_F(LinearAlgebraTests,Lapack_EigenSolverDenseComplex)
{
    EigenTester<dcmplx>(new oml::LapackEigenSolver<dcmplx>(),Neigen).RunTests();
}


//TEST_F(LinearAlgebraTests,Arpack_EigenSolverDenseReal)
//{
//    EigenTester<double>(new ArpackEigenSolver<double>(),Neigen,Neigen-2).RunTests();
//}
//
//TEST_F(LinearAlgebraTests,Arpack_EigenSolverDenseComplex)
//{
//    EigenTester<dcmplx>(new ArpackEigenSolver<dcmplx>(),Neigen,Neigen-2).RunTests();
//}
//
//
//TEST_F(LinearAlgebraTests,Arpack_EigenSolverSparseReal)
//{
//    SparseEigenTester<double>(new ArpackEigenSolver<double>(),Neigen,Neigen-2).RunTests();
//}
//
//TEST_F(LinearAlgebraTests,Arpack_EigenSolverSparseComplex)
//{
//    SparseEigenTester<dcmplx>(new ArpackEigenSolver<dcmplx>(),Neigen,Neigen-2).RunTests();
//}
//
//

//
//TEST_F(LinearAlgebraTests,oml_SVDRandomComplexMatrix_10x10)
//{
//    SetupC(10,10);
//    auto [U,s,Vdagger]=oml_CSVDecomp(itsAC); //Solve A=U*s*conj(V)
//
//    MatrixCT V=Transpose(conj(Vdagger));
//    EXPECT_NEAR(Max(fabs(Transpose(conj(U))*U-itsIC)),0.0,eps);
//    EXPECT_NEAR(Max(fabs(V*Vdagger-itsIC)),0.0,eps);
//    EXPECT_NEAR(Max(fabs(U*s*Vdagger-itsAC)),0.0,eps);
//}
//
//
//TEST_F(LinearAlgebraTests,oml_SVDRandomComplexMatrix_10x5)
//{
//    SetupC(10,5);
//    auto [U,s,Vdagger]=oml_CSVDecomp(itsAC);
//
//    MatrixCT V=Transpose(conj(Vdagger));
////    EXPECT_NEAR(Max(fabs(Transpose(conj(U))*U-itsIC)),0.0,eps); //Not true for 10x5
//    EXPECT_NEAR(Max(fabs(V*Vdagger-itsIC)),0.0,eps);
//    EXPECT_NEAR(Max(fabs(U*s*Vdagger-itsAC)),0.0,eps);
//}
//
//
//TEST_F(LinearAlgebraTests,oml_SVDRandomComplexMatrix_5x10)
//{
//    SetupC(5,10);
//    auto [U,s,Vdagger]=oml_CSVDecomp(itsAC);
//
//    MatrixCT V=Transpose(conj(Vdagger));
//    EXPECT_NEAR(Max(fabs(Transpose(conj(U))*U-itsIC)),0.0,eps);
//    //EXPECT_NEAR(Max(fabs(V*Vdagger-itsIC)),0.0,eps); Not true for 5x10
//    EXPECT_NEAR(Max(fabs(U*s*Vdagger-itsAC)),0.0,eps);
//}
//
//TEST_F(LinearAlgebraTests,OML_EigenSolverComplexHermitian_oldUI)
//{
//    SetupH(50);
//    MatrixCT Acopy(itsAC);
//    int ierr=0;
//    ch(itsAC, itsWR ,true,ierr);
//    EXPECT_EQ(ierr,0);
//    MatrixCT diag=Transpose(conj(itsAC))*Acopy*itsAC;
//    for (int i=1;i<=itsWR.size();i++) diag(i,i)-=itsWR(i);
//    EXPECT_NEAR(Max(fabs(diag)),0.0,Neigen*eps);
//}
//
//TEST_F(LinearAlgebraTests,OML_EigenSolverComplexHermitian)
//{
//    SetupH(50);
//    auto [U,w]=oml_Diagonalize(itsAC);
//    MatrixCT diag=Transpose(conj(U))*itsAC*U;
//    for (int i=1;i<=w.size();i++) diag(i,i)-=w(i);
//    EXPECT_NEAR(Max(fabs(diag)),0.0,Neigen*eps);
//}
//
//TEST_F(LinearAlgebraTests,omlDiagonalMatrix_double)
//{
//    int N=10;
//    Vector<double> v(N);
//    Fill(v,-1.);
//    DiagonalMatrixRT d(v);
//    MatrixRT M(N,N);
//    FillRandom(M);
//    MatrixRT Md=M*d;
//    EXPECT_NEAR(Max(fabs(M+Md)),0.0,eps);
//    MatrixRT dM=d*M;
//    EXPECT_NEAR(Max(fabs(M+dM)),0.0,eps);
//    MatrixRT dMdMdM=d*M*d*M*d*M;
//    EXPECT_NEAR(Max(fabs(M*M*M+dMdMdM)),0.0,eps);
//
//}
//
//TEST_F(LinearAlgebraTests,omlDiagonalMatrix_complex)
//{
//    int N=10;
//    Vector<dcmplx> v(N);
//    Fill(v,dcmplx(-1.0));
//    DiagonalMatrix<dcmplx> d(v);
//    MatrixCT M(N,N);
//    FillRandom(M);
//    MatrixCT Md=M*d;
//    EXPECT_NEAR(Max(fabs(M+Md)),0.0,eps);
//    MatrixCT dM=d*M;
//    EXPECT_NEAR(Max(fabs(M+dM)),0.0,eps);
//    MatrixCT dMdMdM=d*M*d*M*d*M;
//    EXPECT_NEAR(Max(fabs(M*M*M+dMdMdM)),0.0,eps);
//
//}
//TEST_F(LinearAlgebraTests,omlDiagonalMatrix_complex_double)
//{
//    int N=10;
//    Vector<double> v(N);
//    Fill(v,-1.0);
//    DiagonalMatrix<double> d(v);
//    MatrixCT M(N,N);
//    FillRandom(M);
//    MatrixCT Md=M*d;
//    EXPECT_NEAR(Max(fabs(M+Md)),0.0,eps);
//    MatrixCT dM=d*M;
//    EXPECT_NEAR(Max(fabs(M+dM)),0.0,eps);
//    MatrixCT dMdMdM=d*M*d*M*d*M;
//    EXPECT_NEAR(Max(fabs(M*M*M+dMdMdM)),0.0,eps);
//
//}
//
//TEST_F(LinearAlgebraTests,LapackQRSolverReal)
//{
//    QRTester<double>(new LapackQRSolver<double>(),Nqr  ,Nqr  ).RunTests();
//    QRTester<double>(new LapackQRSolver<double>(),Nqr/2,Nqr  ).RunTests();
//    QRTester<double>(new LapackQRSolver<double>(),Nqr  ,Nqr/2).RunTests();
//}
//
//TEST_F(LinearAlgebraTests,LapackQRSolverComplex)
//{
//    QRTester<dcmplx>(new LapackQRSolver<dcmplx>(),Nqr  ,Nqr  ).RunTests();
//    QRTester<dcmplx>(new LapackQRSolver<dcmplx>(),Nqr/2,Nqr  ).RunTests();
//    QRTester<dcmplx>(new LapackQRSolver<dcmplx>(),Nqr  ,Nqr/2).RunTests();
//}
//
//TEST_F(LinearAlgebraTests,LapackRankRevealingQR_Real)
//{
//    QRTester<double>(new LapackQRSolver<double>(),Nqr  ,Nqr  ).RunRRTests(Nqr/2);
//    QRTester<double>(new LapackQRSolver<double>(),Nqr-1,Nqr  ).RunRRTests(Nqr/2);
//    QRTester<double>(new LapackQRSolver<double>(),Nqr  ,Nqr-1).RunRRTests(Nqr/2);
//}
//
//TEST_F(LinearAlgebraTests,LapackRankRevealingQR_Complex)
//{
//    QRTester<dcmplx>(new LapackQRSolver<dcmplx>(),Nqr  ,Nqr  ).RunRRTests(Nqr/2);
//    QRTester<dcmplx>(new LapackQRSolver<dcmplx>(),Nqr-1,Nqr  ).RunRRTests(Nqr/2);
//    QRTester<dcmplx>(new LapackQRSolver<dcmplx>(),Nqr  ,Nqr-1).RunRRTests(Nqr/2);
//}
//
//
//TEST_F(LinearAlgebraTests,LapackLinearSolverReal)
//{
//    LinearTester<double>(new LapackLinearSolver<double>(),Nlinear).RunTests();
//}
//
//TEST_F(LinearAlgebraTests,LapackLinearSolverComplex)
//{
//    LinearTester<dcmplx>(new LapackLinearSolver<dcmplx>(),Nlinear).RunTests();
//}
//

TEST_F(LinearAlgebraTests,Lapack_Cholsky)
{
    CholskyTester(Ncholsky).RunTests();
}

