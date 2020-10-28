// File: UT-Matrix-double.cc  Unit test the Matrix class for double data types.

// Copyright (1994-2020), Jan N. Reimers

//#define WARN_DEEP_COPY 1
#include "oml/matrix.h"
#include "oml/vector.h"
#include "oml/random.h"
#include "stopw.h"
#include "gtest/gtest.h"
#include <omp.h>
#include <iostream>
#include <fstream>
#include <complex>
#include <vector>

using std::cout;
using std::endl;

template <class T> class BenchmarkTests : public ::testing::Test
{
public:
    BenchmarkTests()
    {
        StreamableObject::SetToPretty();
    }
};

class BenchmarkRealTests : public ::testing::Test
{
public:
    BenchmarkRealTests()
    {
        StreamableObject::SetToPretty();
    }
};

Matrix<double> Ag;

inline Matrix<double> NullFunction(int i, int N)
{
    return Ag; //move should be automatic
}
inline Matrix<double> NullFunction_move(int i, int N)
{
    return std::move(Ag); //move should be automatic
}

std::vector<double> Vg;
inline std::vector<double>  NullFunctionV(int i, int N)
{
    return Vg;
}
inline std::vector<double>  NullFunctionVmove(int i, int N)
{
    return std::move(Vg);
}


TEST_F(BenchmarkRealTests,MatrixCopy)
{
    index_t N=400,Nreplicates=10000;

    Matrix<double> A(N,N),B(N,N),C(N,N);
    Ag.SetLimits(N,N);
    Vg.resize(N*N);
    FillRandom(A);
    FillRandom(B);
    Fill(C,0.0);
    #ifdef DEBUG
        double expected_matrix_shallow_copy_time=0.06   ;
        double expected_matrix_deep_copy_time=5.0;
        double expected_shallow_copy_time=0.0001;
        double expected_deep_copy_time=4.0;
    #else
        double expected_matrix_shallow_copy_time=0.0007   ;
        double expected_matrix_deep_copy_time=5.0;
        double expected_shallow_copy_time=0.0004;
        double expected_deep_copy_time=4.0;
    #endif
    double matrix_shallow_copy_time;
    double matrix_deep_copy_time;
    double vector_shallow_copy_time;
    double vector_deep_copy_time;


    cout << "------------------ Force Deep copy ------------------------" << endl;
    double t_start,t_stop,tsum=0;
    for (int i=1;i<=Nreplicates;i++)
    {
        t_start=omp_get_wtime();
        B=NullFunction(i,N);
        Ag(1,1)=5.0; //Force copy on write == deep copy
        t_stop=omp_get_wtime();
        tsum+=t_stop-t_start;
    }
    matrix_deep_copy_time=tsum/Nreplicates*1000/(N*N)*1000000;
    cout << "Matrix deep copy time=" << matrix_deep_copy_time << "ms/Mbyte" << endl;
    EXPECT_LT(matrix_deep_copy_time,expected_matrix_deep_copy_time);


    tsum=0.0;
    {
        cout << "------------------ Matrix shallow copy ------------------------" << endl;
        tsum=0.0;
        for (int i=1;i<=Nreplicates;i++)
        {
            Ag.SetLimits(N,N);
            t_start=omp_get_wtime();
     //       Matrix<double> B1(A);
//            Matrix<double> B1=NullFunction(i,N);
            B=NullFunction_move(i,N) ;
            //B1(1,1)=5.0; //Force copy on write == deep copy
            t_stop=omp_get_wtime();
            tsum+=t_stop-t_start;
        }
        matrix_shallow_copy_time=tsum/Nreplicates*1000/(N*N)*1000000;
        cout << "Matrix copy time=" << matrix_shallow_copy_time << "ms/Mbyte" << endl;
        EXPECT_LT(matrix_shallow_copy_time,expected_matrix_shallow_copy_time);
    }
    {
        cout << "------------------std::vector deep copy ------------------------" << endl;
        tsum=0.0;
        std::vector<double> v(N*N);
        for (int i=1;i<=Nreplicates;i++)
        {
            Vg.resize(N*N);
            t_start=omp_get_wtime();
//            std::vector<double> B1=NullFunctionV(i,N); this is shallow too!
            v=NullFunctionV(i,N);
            t_stop=omp_get_wtime();
            tsum+=t_stop-t_start;
        }
        vector_deep_copy_time=tsum/Nreplicates*1000/(N*N)*1000000;
        cout << "vector deep copy time=" << vector_deep_copy_time << "ms/Mbyte" << endl;
        EXPECT_LT(vector_deep_copy_time,expected_deep_copy_time);
    }
     EXPECT_LT(matrix_shallow_copy_time/matrix_deep_copy_time,0.001);
     EXPECT_LT(vector_shallow_copy_time/vector_deep_copy_time,0.001);
    {
        cout << "------------------std::vector shallow copy ------------------------" << endl;
        tsum=0.0;
        for (int i=1;i<=Nreplicates;i++)
        {
            Vg.resize(N*N);
            t_start=omp_get_wtime();
            std::vector<double> B1=NullFunctionVmove(i,N);
            t_stop=omp_get_wtime();
            tsum+=t_stop-t_start;
        }
        vector_shallow_copy_time=tsum/Nreplicates*1000/(N*N)*1000000;
        cout << "vector shallow copy time=" << vector_shallow_copy_time << "ms/Mbyte" << endl;
        EXPECT_LT(vector_shallow_copy_time,expected_shallow_copy_time);
    }
}



template <class T> void mmul(     Matrix<T>& C,const Matrix<T>& A, const Matrix<T>& B);
template <class T> T    vmul(const Vector<T> V1,const Matrix<T>& A, const Vector<T>& V2);

template <class T> void mmul_nomp(Matrix<T>& C,const Matrix<T>& A, const Matrix<T>& B)
{
    typename Matrix<T>::Subscriptor s(C);
    for (index_t i : C.rows())
        for (index_t j : C.cols())
        {
            double t=0.0;
            for (index_t k : A.cols())
                t+=A(i,k)*B(k,j);
            s(i,j)=t;
        }
}

template <class T> void mmul_mp(Matrix<T>& C,const Matrix<T>& A, const Matrix<T>& B)
{
    {
        typename Matrix<T>::Subscriptor s(C);
        index_t N=C.GetNumRows();
//  OpenMP can't parse this range loop yet.
//        for (index_t i : A.rows())
        #pragma omp parallel for collapse(2)
          for (index_t i=A.GetLimits().Row.Low;i<=A.GetLimits().Row.High;i++)
            for (index_t j=B.GetLimits().Col.Low;j<=B.GetLimits().Col.High;j++)
            {
                double t=0.0;
                for (index_t k=1;k<=N;k++)
                    t+=A(i,k)*B(k,j);
                s(i,j)=t;
            }
    }
}

TEST_F(BenchmarkRealTests,OpenMPParallel)
{
#ifdef DEBUG
    index_t N=20,Nreplicates=10;
#else
    index_t N=500,Nreplicates=10;
#endif

    Matrix<double> A(N,N),B(N,N),C(N,N);
    FillRandom(A);
    FillRandom(B);
    Fill(C,0.0);

    //omp_set_num_threads(8);

    #pragma omp parallel
    {
        #pragma omp single
        std::cout << "Number of available threads: " << omp_get_num_threads() << std::endl;
    }
    double MFlops_nomp(0.0),MFlops_mp(0.0),MFlops_exp(0.0),t_start,t_stop;
    for (int i=1;i<=Nreplicates;i++)
    {
        t_start=omp_get_wtime();
        mmul_nomp(C,A,B);
        t_stop=omp_get_wtime();
        MFlops_nomp+=1e-6*N*N*N*2/(t_stop-t_start);

        t_start=omp_get_wtime();
        mmul_mp(C,A,B);
        t_stop=omp_get_wtime();
        MFlops_mp+=1e-6*N*N*N*2/(t_stop-t_start);

        t_start=omp_get_wtime();
        C=A*B;
        t_stop=omp_get_wtime();
        MFlops_exp+=1e-6*N*N*N*2/(t_stop-t_start);
    }
    MFlops_nomp/=Nreplicates;MFlops_mp/=Nreplicates ;MFlops_exp/=Nreplicates;
    std::cout << "Hand coded   No   MP " << MFlops_nomp << " MFlops" << std::endl;
    std::cout << "Hand coded   With MP " << MFlops_mp << " MFlops" << std::endl;
    std::cout << "Exp Template With MP " << MFlops_exp << " MFlops" << std::endl;
    double mp_improvement=MFlops_mp/MFlops_nomp ;
    std::cout << "Ratio   " << mp_improvement << std::endl;
    EXPECT_GT(mp_improvement,3.0);
    EXPECT_NEAR(MFlops_mp/MFlops_exp,1.0,.05);
}

TYPED_TEST_SUITE_P(BenchmarkTests);

TYPED_TEST_P(BenchmarkTests,MatrixAlgebraPerformance)
{
    typedef Matrix<TypeParam> MatrixT;
    typedef Vector <TypeParam> VectorT;
#ifdef DEBUG
    index_t N=20,Nreplicates=10;
#else
    int N=200,Nreplicates=20;
#endif
    MatrixT A(N,N),B(N,N),C(N,N);
    VectorT V(N);
    FillRandom<TypeParam>(A);
    FillRandom<TypeParam>(B);
    FillRandom<TypeParam>(V);
    double t_start=0,t_stop=0;
    std::vector<TypeParam> ds(Nreplicates+1);
    int Ntype=sizeof(TypeParam)/8;
    cout.precision(1);

    double Flops=0.0;
    for (int i=1;i<=Nreplicates;i++)
    {
        t_start=omp_get_wtime();
        C=A*B;
        t_stop=omp_get_wtime();
        Flops+=Ntype*N*N*N*2/(t_stop-t_start);
    }
    std::cout << "A*B   Exp  template " << std::setw(12) << std::fixed << Flops*1e-6/Nreplicates << " MFlops" << std::endl;

    Flops=0.0;
    for (int i=1;i<=Nreplicates;i++)
    {
        t_start=omp_get_wtime();
        mmul(C,A,B);
        t_stop=omp_get_wtime();
        Flops+=Ntype*N*N*N*2/(t_stop-t_start);
    }
    std::cout << "A*B   Hand coaded   " << std::setw(12) << Flops*1e-6/Nreplicates << " MFlops" << std::endl;


    N=1000;
    Nreplicates=50;
    A.SetLimits(N,N);
    V.SetLimits(N);
    FillRandom<TypeParam>(A);
    FillRandom<TypeParam>(V);
    ds.resize(Nreplicates+1);

    Flops=0.0;
    for (int i=1;i<=Nreplicates;i++)
    {
        t_start=omp_get_wtime();
        ds[i]=V*A*V;
        t_stop=omp_get_wtime();
        Flops+=Ntype*2.*(N*N+N)/(t_stop-t_start);
    }
    std::cout << "V*A*V Exp  template " << std::setw(12) << Flops*1e-6/Nreplicates << " MFlops" << std::endl;
    Flops=0.0;
    for (int i=1;i<=Nreplicates;i++)
    {
        t_start=omp_get_wtime();
        ds[i]=vmul(V,A,V);
        t_stop=omp_get_wtime();
        Flops+=Ntype*2.*(N*N+N)/(t_stop-t_start);
    }
    std::cout << "V*A*V Hand coaded 1 " << std::setw(12) << Flops*1e-6/Nreplicates << " MFlops" << std::endl;

    Flops=0.0;
    for (int i=1;i<=Nreplicates;i++)
    {
        t_start=omp_get_wtime();
        ds[i]=vmul1(V,A,V);
        t_stop=omp_get_wtime();
        Flops+=Ntype*2.*(N*N+N)/(t_stop-t_start);
    }
    std::cout << "V*A*V Hand coaded 2 " << std::setw(12) << Flops*1e-6/Nreplicates << " MFlops" << std::endl;

}

template <class T> void mmul(Matrix<T>& C,const Matrix<T>& A, const Matrix<T>& B)
{
  int N=A.GetLimits().Row.High;
  typename Matrix<T>::Subscriptor s(C);
  for (int i=1;i<=N;i++)
    for (int j=1;j<=N;j++)
    {
      T t(0);
      for (int k=1;k<=N;k++) t+=A(i,k)*B(k,j);
      s(i,j)=t;
    }
}

template <class T> T vmul(const Vector<T> V1,const Matrix<T>& A, const Vector<T>& V2)
{
  T ret=0;
  int N=A.GetLimits().Row.High;
  for (int i=1;i<=N;i++)
  {
    T t(0);
    for (int j=1;j<=N;j++)
        t+=A(i,j)*V2(j);
    ret+=V1(i)*t;
  }
  return ret;
}
template <class T> T vmul1(const Vector<T> V1,const Matrix<T>& A, const Vector<T>& V2)
{
  T ret=0;
  int N=A.GetLimits().Row.High;
  for (int j=1;j<=N;j++)
  {
    T t(0);
    for (int i=1;i<=N;i++)
        t+=V1(i)*A(i,j);
    ret+=V2(j)*t;
  }
  return ret;
}



REGISTER_TYPED_TEST_SUITE_P(BenchmarkTests,
            MatrixAlgebraPerformance
            );
//using MyTypes = ::testing::Types<double>;
using MyTypes = ::testing::Types<double,std::complex<double>>;
INSTANTIATE_TYPED_TEST_SUITE_P(My, BenchmarkTests, MyTypes);

