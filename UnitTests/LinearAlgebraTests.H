#ifndef LINEARALGEBRATESTS_H_INCLUDED
#define LINEARALGEBRATESTS_H_INCLUDED

//#include "NumericalMethods/SparseSVDSolver.H"
#include "oml/numeric/EigenSolver.H"
#include "oml/numeric/SVDSolver.H"
//#include "NumericalMethods/SparseEigenSolver.H"
#include "oml/numeric/LapackEigenSolver.H"
#include "oml/numeric/LapackSVDSolver.H"
//#include "NumericalMethods/QRSolver.H"
//#include "NumericalMethods/LinearSolver.H"
#include "oml/random.h"
#include "oml/matrix.h"
#include "oml/diagonalmatrix.h"
#include "oml/fakedouble.h"
#include <iostream>
#include <complex>

using std::cout;
using std::endl;

typedef std::complex<double> dcmplx;

//
//  Randomly put zeros into a dense matrix until a request density is achieved
//
template <class T> void MakeSparse(Matrix<T>& m,double TargetDensity)
{
    int Ntotal=m.size(),N=m.size();
    int Nr=m.GetNumRows(),Nc=m.GetNumCols();
    double density=1.0;
    while (density>TargetDensity)
    {
        int ir=static_cast<int>(OMLRand<float>()*Nr)+1;
        int ic=static_cast<int>(OMLRand<float>()*Nc)+1;
        if (m(ir,ic)!=0.0)
        {
            m(ir,ic)=0.0;
            N--;
            density=static_cast<double>(N)/Ntotal;
        }

    }

}

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
        MakeSparse(ret,Density);
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

template <typename T>
class EigenTester : public EigenBase<T>
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


//template <typename T>
//class SparseEigenTester
//    : public EigenBase<T>
//    , public virtual SparseEigenSolverClient<T>
//{
//    using Base=EigenBase<T>;
//    using Base::N;
//    using Base::Nev;
//    using Base::eps;
//    using Base::I;
//public:
//    SparseEigenTester(SparseEigenSolver<T>* s, int n,double density=1.0)
//        : Base(n,density)
//        , itsSolver(s)
//        , A(n,n)
//    {}
//    SparseEigenTester(SparseEigenSolver<T>* s, int n,int Nev,double density=1.0)
//        : Base(n,Nev,density)
//        , itsSolver(s)
//        , A(n,n)
//    {}
//    ~SparseEigenTester()
//    {
//        delete itsSolver;
//    }
//
//    virtual void RunTests(bool useNormal=false)
//    {
//        A=Base::MakeSymMatrix();
//        {
//            auto [U,d]=itsSolver->Solve(A,eps,Nev);
//            ValidateSymmetric(U,d);
//            itsSolver->Reset();
//        }
//        {
//            auto [U,d]=itsSolver->Solve(this,eps,Nev);
//            ValidateSymmetric(U,d);
//            itsSolver->Reset();
//        }
//        if (Nev/2>=1)
//        {
//            auto [U,d]=itsSolver->Solve(A,eps,Nev/2); //Get only half the eigen values
//            ValidateSymmetric(U,d);
//            itsSolver->Reset();
//        }
//        if (Nev/2>=1)
//        {
//            auto [U,d]=itsSolver->Solve(this,eps,Nev/2); //Get only half the eigen values
//            ValidateSymmetric(U,d);
//            itsSolver->Reset();
//        }
//        if (useNormal)
//            A=Base::MakeNormalMatrix(); //Make a non sym normal matrix
//        else
//            A=Base::MakeNonSymMatrix(); //Make an unrestricted non sym matrix
//        {
//            auto [U,d]=itsSolver->SolveRightNonSym(A,eps,Nev);
//            ValidateRight(U,d);
//            itsSolver->Reset();
//        }
//        {
//            auto [U,d]=itsSolver->SolveLeft_NonSym(A,eps,Nev);
//            ValidateLeft(U,d);
//            itsSolver->Reset();
//        }
//        {
//            auto [U,d]=itsSolver->SolveRightNonSym(this,eps,Nev);
//            ValidateRight(U,d);
//            itsSolver->Reset();
//        }
//        if (Nev/2>=1)
//        {
//            auto [U,d]=itsSolver->SolveRightNonSym(A,eps,Nev/2);
//            ValidateRight(U,d);
//            itsSolver->Reset();
//        }
//        if (Nev/2>=1)
//        {
//            auto [U,d]=itsSolver->SolveLeft_NonSym(A,eps,Nev/2);
//            ValidateLeft(U,d);
//            itsSolver->Reset();
//        }
//        if (Nev/2>=1)
//        {
//            auto [U,d]=itsSolver->SolveRightNonSym(this,eps,Nev/2);
//            ValidateRight(U,d);
//            itsSolver->Reset();
//        }
//    }
//    virtual void ValidateSymmetric(const Matrix<T>& U, const Vector<double>& d)
//    {
//        Matrix<T> UT=Transpose(conj(U));
//        int mn=U.GetNumCols();
//        Matrix<T> Iu(mn,mn);
//        Unit(Iu);
//        //cout << "U*UT=" << UT*U << endl;
//        EXPECT_NEAR(Max(fabs(UT*U-Iu)),0.0,2*N*eps);
//        Matrix<T> diag=UT*A*U;
////        cout << "diag=" << diag << endl;
//        int Ne=d.size();
//        for (int i=1; i<=Ne; i++) diag(i,i)-=d(i);
//        EXPECT_NEAR(Max(fabs(diag)),0.0,2*N*eps);
//    }
//    virtual void ValidateRight(const Matrix<dcmplx>& U, const Vector<dcmplx>& d)
//    {
//        int Ne=d.size();
//        if (Ne>0)
//        {
//            for (int i=1; i<=Ne; i++)
//            {
//                Vector<dcmplx> Ui=U.GetColumn(i);
//                Vector<dcmplx> residuals=A*Ui-d(i)*Ui; //Anonymous warning from here
//                double res=Max(fabs(residuals));
//                EXPECT_NEAR(res,0.0,15*eps);
//            }
//        }
//    }
//    virtual void ValidateLeft(const Matrix<dcmplx>& U, const Vector<dcmplx>& d)
//    {
//        int Ne=d.size();
//        if (Ne>0)
//        {
//            for (int i=1; i<=Ne; i++)
//            {
//                Vector<dcmplx> Ui=U.GetColumn(i);
//                Vector<dcmplx> residuals=Ui*A-d(i)*Ui;
//                double res=Max(fabs(residuals));
//                EXPECT_NEAR(res,0.0,15*eps);
//            }
//        }
//    }
//
//    virtual long int GetSize  () const
//    {
//        return N;
//    }; //2^L Hilbert space size
//    virtual void     DoMatVecContraction (int N, const T* xvec, T* yvec) const
//    {
//        A.DoMVMultiplication(N,xvec,yvec);
//    };
//
//
//    SparseEigenSolver<T>* itsSolver;
//    SparseMatrix     <T>  A;
//};

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
        MakeSparse(ret,Density);
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


template <typename T>
class SVDTester : public SVDBase<T>
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

#include "oml/smatrix.h"
#include "oml/numeric/LapackCholsky.H"
#include "oml/numeric.h"

class CholskyTester 
{
public:
    CholskyTester(int m)
        : M(m), eps(1e-13)
        , A(m)
    {}

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
    
    virtual void RunTests()
    {
        A=MakeMatrix();
//        cout << "A=" << A << endl;
        {
            auto V=oml::LapackCholsky(A);
//            cout << "V=" << V << endl;
            Validate(V,A);
            
            Matrix<double> Voml(A);
            Cholsky(Voml);
//            cout << "Voml=" << Voml << endl;
            Validate(Transpose(Voml),A);
        }
    }
    virtual void Validate(const Matrix<double>& V, const SMatrix<double>& A1)
    {
//        cout << Matrix<double>(V*Transpose(V)) << endl;
//        cout << Matrix<double>(Transpose(V)*V) << endl;
        Matrix<double> dA=A1-Transpose(V)*V;
        EXPECT_NEAR(Max(fabs(dA)),0.0,2*M*eps);
    }

    int M;
    double eps;
    SMatrix<double>  A;
};



//template <typename T>
//class SparseSVDTester
//    : public SVDBase<T>
//    , public virtual SparseSVDSolverClient<T>
//{
//    using Base=SVDBase<T>;
//    using DMatRT=typename Base::DMatRT;
//    using Base::M;
//    using Base::N;
//    using Base::eps;
//    using Base::I;
//
//public:
//    SparseSVDTester(SparseSVDSolver<T>* s, int m, int n,double density=1.0)
//        : Base(m,n,density)
//        , itsSolver(s)
//        , A(m,n)
//    {}
//    ~SparseSVDTester()
//    {
//        delete itsSolver;
//    }
//
//    virtual void RunTests()
//    {
//        A=Base::MakeMatrix();
//        {
//            auto [U,s,VT]=itsSolver->Solve(A,eps,Min(M,N));
//            Validate(U,s,VT);
//        }
//        {
//            //Client version
//            auto [U,s,VT]=itsSolver->Solve(this,eps,Min(M,N));
//            Validate(U,s,VT);
//        }
//    }
//    virtual void Validate(const Matrix<T>& U, const DMatRT& s,const Matrix<T>& VT)
//    {
//        Base::Validate(U,s,VT); //Check U,V orthonormality
//        Matrix<T> A1=U*s*VT;
//        Matrix<T> dA=A-A1;
//        EXPECT_NEAR(Max(fabs(dA)),0.0,2*sqrt(N*M)*eps);
//    }
////
//// SparseSVDSolverClient interface
////
//    virtual int  GetNumRows() const
//    {
//        return M;
//    }
//    virtual int  GetNumCols() const
//    {
//        return N;
//    }
//    virtual void DoMatVecContraction (int N, const T* xvec, T* yvec,int transpose) const
//    {
//        A.DoMVMultiplication(M,N,xvec,yvec,transpose);
//    }
//
//private:
//    SparseSVDSolver<T>* itsSolver;
//    SparseMatrix   <T>  A;
//};
//
//
////
////   QR Tester classes
////
//template <typename T> class QRBase
//{
//public:
//
//    QRBase(int m, int n,double density=1.0)
//        : M(m), N(n), Density(density),eps(1e-13)
//        , I(Min(M,N),Min(M,N))
//    {
//        Unit(I);
//    }
//    Matrix<T> MakeMatrix()
//    {
//        Matrix<T> ret(M,N);
//        FillRandom(ret);
//        MakeSparse(ret,Density);
//        return ret;
//    }
//    Matrix<T> MakeMatrix(int Rank)
//    {
//        assert(Rank<=M);
//        assert(Rank<=N);
//        Matrix<T> ret(M,N);
//        Fill(ret,T(0.0));
//        for (int i=1;i<=Rank;i++)
//        {
//            Vector<T> v1(M),v2(N);
//            FillRandom(v1);
//            FillRandom(v2);
//            Matrix<T> op=OuterProduct(v1,v2);
//            ret+=op;
//        }
//        MakeSparse(ret,Density);
//        return ret;
//    }
//
//    virtual void Validate(const Matrix<T>& Q)
//    {
//
//    }
//
//
//protected:
//    int M,N;
//    double Density;
//    double eps;
//    Matrix<T> I;
//};
//
//template <typename T>
//class QRTester : public QRBase<T>
//{
//    using Base=QRBase<T>;
//    using Base::M;
//    using Base::N;
//    using Base::eps;
//    using Base::I;
//public:
//    QRTester(QRSolver<T>* s, int m, int n,double density=1.0)
//        : Base(m,n,density)
//        , itsSolver(s)
//        , A(m,n)
//    {}
//    ~QRTester()
//    {
//        delete itsSolver;
//    }
//
//    virtual void RunTests()
//    {
//        A=Base::MakeMatrix();
//        {
//            auto [Q,R]=itsSolver->SolveThinQR(A);
//            EXPECT_TRUE(IsUpperTriangular(R,eps));
//            ValidateQR(Q,R);
//        }
//        {
//            auto [R,Q]=itsSolver->SolveThinRQ(A);
//            EXPECT_TRUE(IsUpperTriangular(R,eps));
//            ValidateRQ(R,Q);
//        }
//
//        {
//            auto [Q,L]=itsSolver->SolveThinQL(A);
//            EXPECT_TRUE(IsLowerTriangular(L,eps));
//            ValidateQR(Q,L);
//        }
//        {
//            auto [L,Q]=itsSolver->SolveThinLQ(A);
//            EXPECT_TRUE(IsLowerTriangular(L,eps));
//            ValidateRQ(L,Q);
//        }
//
//    }
//    //
//    //  Rank revealing QR tests
//    //
//    virtual void RunRRTests(int Rank)
//    {
//        A=Base::MakeMatrix(Rank);
//        {
//            auto [Q,R]=itsSolver->SolveRankRevealingQR(A,10*sqrt(N*M)*eps); //Need a little lattitude to remove all rows.
//            EXPECT_TRUE(IsUpperTriangular(R,eps));
//            ValidateQR(Q,R,Rank);
//        }
//        {
//            auto [R,Q]=itsSolver->SolveRankRevealingRQ(A,10*sqrt(N*M)*eps); //Need a little lattitude to remove all rows.
//            EXPECT_TRUE(IsUpperTriangular(R,eps));
//            ValidateRQ(R,Q,Rank);
//        }
//        {
//            auto [Q,L]=itsSolver->SolveRankRevealingQL(A,10*sqrt(N*M)*eps); //Need a little lattitude to remove all rows.
//            EXPECT_TRUE(IsLowerTriangular(L,eps));
//            ValidateQR(Q,L,Rank);
//        }
//        {
//            auto [L,Q]=itsSolver->SolveRankRevealingLQ(A,10*sqrt(N*M)*eps); //Need a little lattitude to remove all rows.
//            EXPECT_TRUE(IsLowerTriangular(L,eps));
//            ValidateRQ(L,Q,Rank);
//        }
//    }
//
//    virtual void ValidateQR(const Matrix<T>& Q, const Matrix<T>& R)
//    {
//
//        Matrix<T> QT=Transpose(conj(Q));
//        {
//            int mn=Q.GetNumCols();
//            Matrix<T> Iq(mn,mn);
//            Unit(Iq);
//            EXPECT_NEAR(Max(fabs(QT*Q-Iq)),0.0,2*sqrt(N*M)*eps);
//        }
//        if (N==M && Q.GetNumCols()==M)
//        {
//            int mn=Q.GetNumRows();
//            Matrix<T> Iq(mn,mn);
//            Unit(Iq);
//            EXPECT_NEAR(Max(fabs(Q*QT-Iq)),0.0,2*sqrt(N*M)*eps);
//        }
//        Matrix<T> A1=Q*R;
//        Matrix<T> dA=A-A1;
//        EXPECT_NEAR(Max(fabs(dA)),0.0,2*sqrt(N*M)*eps);
//    }
//    virtual void ValidateRQ(const Matrix<T>& R, const Matrix<T>& Q)
//    {
//
//        Matrix<T> QT=Transpose(conj(Q));
//        {
//            int mn=Q.GetNumRows();
//            Matrix<T> Iq(mn,mn);
//            Unit(Iq);
//            EXPECT_NEAR(Max(fabs(Q*QT-Iq)),0.0,2*sqrt(N*M)*eps);
//        }
//        if (N==M && Q.GetNumRows()==M)
//        {
//            int mn=Q.GetNumCols();
//            Matrix<T> Iq(mn,mn);
//            Unit(Iq);
//            EXPECT_NEAR(Max(fabs(QT*Q-Iq)),0.0,2*sqrt(N*M)*eps);
//        }
//        Matrix<T> A1=R*Q;
//        Matrix<T> dA=A-A1;
//        EXPECT_NEAR(Max(fabs(dA)),0.0,2*sqrt(N*M)*eps);
//    }
//    virtual void ValidateQR(const Matrix<T>& Q, const Matrix<T>& R,int Rank)
//    {
//        ValidateQR(Q,R);
//        if (M<=30 && N<=30)
//        {
//            EXPECT_EQ(Q.GetNumCols(),Rank);
//            EXPECT_EQ(R.GetNumRows(),Rank);
//        }
//        else
//        {
//            EXPECT_LE(Q.GetNumCols(),Rank+3); //FOr large M,N it is very has to build a test matrix with an exact Rank
//            EXPECT_LE(R.GetNumRows(),Rank+3);
//        }
//    }
//    virtual void ValidateRQ(const Matrix<T>& R, const Matrix<T>& Q,int Rank)
//    {
//        ValidateRQ(R,Q);
//        if (M<=30 && N<=30)
//        {
//            EXPECT_EQ(Q.GetNumRows(),Rank);
//            EXPECT_EQ(R.GetNumCols(),Rank);
//        }
//        else
//        {
//            EXPECT_LE(Q.GetNumRows(),Rank+3);
//            EXPECT_LE(R.GetNumCols(),Rank+3);
//        }
//    }
//
//    QRSolver<T>* itsSolver;
//    Matrix  <T>  A;
//};
//
//template <typename T> class LinearTester
//{
//public:
//    LinearTester(LinearSolver<T>* s, int n)
//        : itsSolver(s)
//        , N(n)
//        , eps(1e-14)
//        , A(n,n)
//        , b(n)
//    {}
//    ~LinearTester()
//    {
//        delete itsSolver;
//    }
//
//    virtual void RunTests()
//    {
//        A=MakeMatrix('U');
//        FillRandom(b);
//        {
//            Vector<T> xr=itsSolver->SolveUpperTri(A,b); //Solve A*x=b
//            Vector<T> xl=itsSolver->SolveUpperTri(b,A); //Solve b=x*A
//            Validate(xr,xl);
//        }
//        A=MakeMatrix('L');
//        {
//            Vector<T> xr=itsSolver->SolveLowerTri(A,b); //Solve A*x=b
//            Vector<T> xl=itsSolver->SolveLowerTri(b,A); //Solve b=x*A
//            Validate(xr,xl);
//        }
//        A=MakeMatrix('F');
//        {
//            Vector<T> xr=itsSolver->Solve(A,b); //Solve A*x=b
//            Vector<T> xl=itsSolver->Solve(b,A); //Solve b=x*A
//            Validate(xr,xl);
//        }
//    }
//    virtual void Validate(const Vector<T>& xr,const Vector<T>& xl)
//    {
//        EXPECT_NEAR(Max(fabs(b-A*xr)),0.0,sqrt(N)*eps);
//        EXPECT_NEAR(Max(fabs(b-xl*A)),0.0,sqrt(N)*eps);
//    }
//
//    Matrix<T> MakeMatrix(char ul)
//    {
//        Matrix<T> ret(N,N);
//        FillRandom(ret);
//        if (ul=='L')
//        {
//            for (index_t i:ret.rows())
//                for (index_t j:ret.cols(i+1))
//                    ret(i,j)=T(0.0);
//        }
//        else if(ul=='U')
//        {
//            for (index_t j:ret.cols())
//                for (index_t i:ret.rows(j+1))
//                    ret(i,j)=T(0.0);
//        }
//        else if (ul=='F')
//        {}//Noop
//        else
//            assert(false);
//
//        for (index_t i:ret.rows())
//            ret(i,i)=T(sqrt(sqrt(N))); //Ensure some numerical stability.
//
//        return ret;
//    }
//
//    LinearSolver<T>* itsSolver;
//    int       N;
//    double    eps;
//    Matrix<T> A;
//    Vector<T> b;
//};
//

#endif // LINEARALGEBRATESTS_H_INCLUDED
