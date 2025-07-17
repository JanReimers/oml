#include "gtest/gtest.h"
#include <cmath>
#include <complex>
#include "Epsilons.H"

import oml.Matrix;
import oml.SMatrix;
import oml.DiagonalMatrix;
import oml.Vector;
import oml.Solvers;
import oml.Lapack;


typedef std::complex<double> dcmplx;
template <typename T> class EigenBase
{
public:
    EigenBase(int n,double density=1.0)
        : N(n), Nev(n),Density(density),eps(1e-13)
        , I(N,N)
    {
        Unit(I);
    }
    EigenBase(int n,int nev,double density=1.0)
        : N(n), Nev(nev),Density(density),eps(1e-13)
        , I(N,N)
    {
        Unit(I);
    }
    Matrix<T> MakeNonSymMatrix()
    {
        Matrix<T> ret(N,N);
        FillRandom(ret);
        // MakeSparse(ret,Density);
        return ret;
    }
    Matrix<T> MakeSymMatrix()
    {
        Matrix<T> ret=MakeNonSymMatrix();
        ret=Matrix<T>(ret+Transpose(conj(ret)));
        return ret;
    }
    Matrix<T> MakeNormalMatrix();

    virtual void ValidateSymmetric(const Matrix<T>& U)
    {
        Matrix<T> UT=Transpose(conj(U));
        //cout << "U*UT=" << UT*U << endl;
        EXPECT_NEAR(Max(fabs(UT*U-I)),0.0,2*N*eps);
    }


protected:
    int N;
    int Nev;
    double Density;
    double eps;
    Matrix<T> I;
};
template <>  Matrix<double> EigenBase<double>::MakeNormalMatrix()
{
//    return MakeSymMatrix();
    Matrix<double> A=MakeSymMatrix();
    auto [U,d]=oml::LapackEigenSolver<double>().SolveAll(A,eps);
    return U;
//    DiagonalMatrix<dcmplx> e(d.size());
//    e(1)=d(1);
//    for (index_t i=2;i<=d.size();i++)
//    {
//        double im=1.01*OMLRandPos<double>();
//        e(i)  =dcmplx(d(i), im);
//        e(i+1)=dcmplx(d(i+1),-im);
//        i++;
//    }
//    Matrix<dcmplx> Ac=U*e*~U;
//    cout << "A=" << Ac << endl;
//    cout << "A*~A=" << Ac*~Ac << endl;
//    cout << "A*~A-~A*A=" << Ac*~Ac-~Ac*Ac << endl;
//    assert(IsNormal(Ac,1e-13));
//    return Ac;
}
template <>  Matrix<dcmplx> EigenBase<dcmplx>::MakeNormalMatrix()
{
    Matrix<dcmplx> A=MakeSymMatrix();
    auto [U,d]=oml::LapackEigenSolver<dcmplx>().SolveAll(A,eps);
    DiagonalMatrix<dcmplx> e(d.size());
    for (index_t i=1;i<=d.size();i++)
    {
        double im=1.00*OMLRand<double>();
        e(i)  =dcmplx(d(i), im);
    }
    Matrix<dcmplx> Ac=U*e*~U;
    assert(IsNormal(Ac,Ac.GetNumRows()*1e-13));
    return Ac;
}

template <typename T> class EigenTester : public EigenBase<T>
{
    using Base=EigenBase<T>;
    using Base::N;
    using Base::Nev;
    using Base::eps;
    using Base::I;
public:
    EigenTester(oml::EigenSolver<T>* s, int n,double density=1.0)
        : Base(n,density)
        , itsSolver(s)
        , A(n,n)
    {}
    EigenTester(oml::EigenSolver<T>* s, int n,int Nev,double density=1.0)
        : Base(n,Nev,density)
        , itsSolver(s)
        , A(n,n)
    {}
    ~EigenTester()
    {
        delete itsSolver;
    }

    virtual void RunTests(bool useNormal=false)
    {
        A=Base::MakeSymMatrix();
        {
            auto [U,d]=itsSolver->Solve(A,eps,Nev);
            ValidateSymmetric(U,d);
            itsSolver->Reset();
        }
        if (Nev/2>=1)
        {
            auto [U,d]=itsSolver->Solve(A,eps,Nev/2); //Get only half the eigen values
            ValidateSymmetric(U,d);
            itsSolver->Reset();
        }
        if (Nev==N)
        {
            auto [U,d]=itsSolver->SolveAll(A,eps);
            ValidateSymmetric(U,d);
            itsSolver->Reset();
        }
        if (useNormal)
            A=Base::MakeNormalMatrix(); //Make a non sym normal matrix
        else
            A=Base::MakeNonSymMatrix(); //Make an unrestricted non sym matrix
        {
            auto [U,d]=itsSolver->SolveRightNonSym(A,eps,Nev);
            ValidateRight(U,d);
            itsSolver->Reset();
        }
        {
            auto [U,d]=itsSolver->SolveLeft_NonSym(A,eps,Nev);
            ValidateLeft(U,d);
            itsSolver->Reset();
        }
        if (Nev/2>=1)
        {
            auto [U,d]=itsSolver->SolveRightNonSym(A,eps,Nev/2);
            ValidateRight(U,d);
            itsSolver->Reset();
        }
        if (Nev/2>=1)
        {
            auto [U,d]=itsSolver->SolveLeft_NonSym(A,eps,Nev/2);
            ValidateLeft(U,d);
            itsSolver->Reset();
        }
        if (Nev==N)
        {
            auto [U,d]=itsSolver->SolveAllRightNonSym(A,eps);
            ValidateRight(U,d);
            itsSolver->Reset();
        }
    }
    virtual void ValidateSymmetric(const Matrix<T>& U, const Vector<double>& d)
    {
        Matrix<T> UT=Transpose(conj(U));
        int mn=U.GetNumCols();
        Matrix<T> Iu(mn,mn);
        Unit(Iu);
        EXPECT_NEAR(Max(fabs(UT*U-Iu)),0.0,2*N*eps);
        Matrix<T> diag=UT*A*U;
//        cout << "diag=" << diag << endl;
        int Ne=d.size();
        for (int i=1; i<=Ne; i++) diag(i,i)-=d(i);
        EXPECT_NEAR(Max(fabs(diag)),0.0,N*eps);
    }
    virtual void ValidateRight(const Matrix<dcmplx>& U, const Vector<dcmplx>& d)
    {
        int Ne=d.size();
        if (Ne>0)
        {
            for (int i=1; i<=Ne; i++)
            {
                Vector<dcmplx> Ui=U.GetColumn(i);
                Vector<dcmplx> residuals=A*Ui-d(i)*Ui;
                double res=Max(fabs(residuals));
                EXPECT_NEAR(res,0.0,U.GetNumRows()*eps);
            }
        }
    }
    virtual void ValidateLeft(const Matrix<dcmplx>& U, const Vector<dcmplx>& d)
    {
        int Ne=d.size();
        if (Ne>0)
        {
            for (int i=1; i<=Ne; i++)
            {
                Vector<dcmplx> Ui=U.GetColumn(i);
                Vector<dcmplx> residuals=Ui*A-d(i)*Ui;
                double res=Max(fabs(residuals));
                EXPECT_NEAR(res,0.0,U.GetNumRows()*eps);
            }
        }
    }

    oml::EigenSolver<T>* itsSolver;
    Matrix     <T>  A;
};


//
//   SVD Tester classes
//
template <typename T> class SVDBase
{
public:
    using DMatRT=DiagonalMatrix<double>;

    SVDBase(int m, int n,double density=1.0)
        : M(m), N(n), Density(density),eps(1e-13)
        , I(Min(M,N),Min(M,N))
    {
        Unit(I);
    }
    Matrix<T> MakeMatrix()
    {
        Matrix<T> ret(M,N);
        FillRandom(ret);
        // MakeSparse(ret,Density);
        return ret;
    }

    virtual void Validate(const Matrix<T>& U, const DMatRT& s,const Matrix<T>& VT)
    {
        if (N<=M)
        {
            Matrix<T> V=Transpose(conj(VT));
            EXPECT_NEAR(Max(fabs(V*VT-I)),0.0,2*sqrt(N*M)*eps);
        }
        if (N>=M)
        {
            Matrix<T> UT=Transpose(conj(U));
            //cout << "U*UT=" << UT*U << endl;
            EXPECT_NEAR(Max(fabs(UT*U-I)),0.0,2*sqrt(N*M)*eps);
        }
    }


protected:
    int M,N;
    double Density;
    double eps;
    Matrix<T> I;
};
template <typename T> class SVDTester : public SVDBase<T>
{
    using Base=SVDBase<T>;
    using DMatRT=typename Base::DMatRT;
    using Base::M;
    using Base::N;
    using Base::eps;
    using Base::I;
public:
    SVDTester(typename oml::SVDSolver<T>* s, int m, int n,double density=1.0)
        : Base(m,n,density)
        , itsSolver(s)
        , A(m,n)
    {}
    ~SVDTester()
    {
        delete itsSolver;
    }

    virtual void RunTests()
    {
        A=Base::MakeMatrix();
        {
            auto [U,s,VT]=itsSolver->Solve(A,Min(M,N));
            Validate(U,s,VT);
        }
        {
            auto [U,s,VT]=itsSolver->SolveAll(A);
            Validate(U,s,VT);
        }
    }
    virtual void Validate(const Matrix<T>& U, const DMatRT& s,const Matrix<T>& VT)
    {
        Base::Validate(U,s,VT); //Check U,V orthonormality
        Matrix<T> A1=U*s*VT;
        Matrix<T> dA=A-A1;
        EXPECT_NEAR(Max(fabs(dA)),0.0,2*sqrt(N*M)*eps);
    }

    oml::SVDSolver<T>* itsSolver;
    Matrix   <T>  A;
};
class CholskyTester 
{
public:
    CholskyTester(int m)
        : M(m), eps(1e-13)
        , A(m), I(m)
    {
        Unit(I);
    }

    SMatrix<double> MakeMatrix()
    {
        SMatrix<double> m(M);
        for (int i=m.GetLimits().Row.Low; i<=m.GetLimits().Row.High; i++)
        for (int j=i; j<=m.GetLimits().Col.High; j++)
        {
            double del=(i-j)/2.0;
            m(i,j)=exp(-del*del);
        }
        return m;
    }
    
    Matrix<double> MakeUpperTriMatrix()
    {
//        Matrix<double> m(M,M);
//        FillRandom(m);
//        for (int i=m.GetLimits().Row.Low; i<=m.GetLimits().Row.High; i++)
//        for (int j=m.GetLimits().Col.Low; j<i; j++)
//            m(i,j)=0.0;
//        return m;
        return oml::LapackCholsky(MakeMatrix());
    }
    
    virtual void RunTests()
    {
        A=MakeMatrix();
        {
            auto V=oml::LapackCholsky(A);
            Validate(V,A);
        }
        {
            auto Ainv=oml::LapackInvertSymmetric(A);
            Validate(Ainv,A);
        }
        auto T=MakeUpperTriMatrix();
        {
            auto Tinv=oml::LapackInvertTriangular(T);
            Validate(Tinv,T);
        }
    }
    virtual void Validate(const Matrix<double>& V, const SMatrix<double>& A1)
    {
        Matrix<double> dA=A1-Transpose(V)*V;
        EXPECT_NEAR(Max(fabs(dA)),0.0,2*M*eps);
    }
    virtual void Validate(const SMatrix<double>& Ainv, const SMatrix<double>& A1)
    {
        Matrix<double> dA=A1*Ainv-I;
        EXPECT_NEAR(Max(fabs(dA)),0.0,2*M*eps);
        dA=Ainv*A1-I;
        EXPECT_NEAR(Max(fabs(dA)),0.0,2*M*eps);
    }
    virtual void Validate(const Matrix<double>& Tinv, const Matrix<double>& T)
    {
        Matrix<double> dT=T*Tinv-I;
        EXPECT_NEAR(Max(fabs(dT)),0.0,2*M*eps);
        dT=Tinv*T-I;
        EXPECT_NEAR(Max(fabs(dT)),0.0,2*M*eps);
    }

    int M;
    double eps;
    SMatrix<double>  A,I;
};



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

    oml::Epsilons  itsEps;
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

TEST_F(LinearAlgebraTests,Lapack_EigenSolverDenseReal)
{
    EigenTester<double>(new oml::LapackEigenSolver<double>(),Neigen).RunTests();
}
TEST_F(LinearAlgebraTests,Lapack_EigenSolverDenseComplex)
{
    EigenTester<dcmplx>(new oml::LapackEigenSolver<dcmplx>(),Neigen).RunTests();
}

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

