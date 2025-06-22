// File:: MatrixTests.C  Test the Matrix class with subscriptors.

#include "gtest/gtest.h"
#include "oml2/Matrix.H"

class MatrixTests : public ::testing::Test
{
public:
    MatrixTests() = default;
    ~MatrixTests() override = default;
};

template <std::ranges::range Range> void print(Range v)
{
    for (auto element : v) {
        std::cout << element << " ";
    }
    std::cout << "\n";
}
auto print2D=[](auto rng){for(auto r:rng)print(r);};

using il=std::initializer_list<double>;

TEST_F(MatrixTests, FullMatrix)
{
    FullMatrix mat(3, 4);
    mat(1, 2) = 5.0;
    EXPECT_EQ(mat(1, 2), 5.0);
    FullMatrix mat1{{1,2,3},
                    {4,5,6},
                    {7,8,9},
                    {10,11,12}};
    EXPECT_EQ(mat1(1, 2), 6.0);
    EXPECT_EQ(mat1(2, 1), 8.0);
    EXPECT_EQ(mat1(0, 0), 1.0);
    EXPECT_EQ(mat1(0, 1), 2.0);
    EXPECT_EQ(mat1(0, 2), 3.0);
    mat1.print();
    print(mat1.row(1)); // Should print: [4,5,6]
    print(mat1.col(2)); // Should print: [3,6,9,12]
    EXPECT_TRUE((mat1.row(1)==il{4,5,6}));
    EXPECT_TRUE((mat1.col(2)==il{3,6,9,12}));
}
TEST_F(MatrixTests, UpperTriangularMatrix)
{
    UpperTriangularMatrix mat(3, 4);
    mat(1, 2) = 5.0;
    EXPECT_EQ(mat(1, 2), 5.0);
    UpperTriangularMatrix mat1{{1,2,3},
                               {0,5,6},
                               {0,0,9},
                               {0,0,0}};
    EXPECT_EQ(mat1(1, 2), 6.0);
    EXPECT_EQ(mat1(0, 0), 1.0);
    EXPECT_EQ(mat1(0, 1), 2.0);
    EXPECT_EQ(mat1(0, 2), 3.0);
    mat1.print();
    print(mat1.row(1)); // Should print: [0,5,6]
    print(mat1.col(2)); // Should print: [3,6,9]
    EXPECT_TRUE((mat1.row(1)==il{5,6}));
    EXPECT_TRUE((mat1.col(2)==il{3,6,9}));
}
TEST_F(MatrixTests, DiagonalMatrix)
{
    DiagonalMatrix mat(3, 4);
    mat(1, 1) = 5.0;
    EXPECT_EQ(mat(1, 1), 5.0);
    DiagonalMatrix mat1{{1,0,0},
                        {0,5,0},
                        {0,0,9},
                        {0,0,0}};
    EXPECT_EQ(mat1(1, 1), 5.0);
    EXPECT_EQ(mat1(2, 2), 9.0);
    EXPECT_EQ(mat1(0, 0), 1.0);
    mat1.print();
    print(mat1.row(1)); // Should print: [0,5,0]
    print(mat1.col(2)); // Should print: [0,9]
    EXPECT_TRUE((mat1.row(1)==il{5}));
    EXPECT_TRUE((mat1.col(2)==il{9}));
}
TEST_F(MatrixTests, TriDiagonalMatrix)
{
    TriDiagonalMatrix mat(3, 3);
    mat(1, 0) = 5.0;
    mat(1, 1) = 6.0;
    mat(1, 2) = 7.0;
    EXPECT_EQ(mat(1, 0), 5.0);
    EXPECT_EQ(mat(1, 1), 6.0);
    EXPECT_EQ(mat(1, 2), 7.0);
    TriDiagonalMatrix mat1{{1,2,0},
                           {5,6,7},
                           {0,8,9}};
    EXPECT_EQ(mat1(1, 0), 5.0);
    EXPECT_EQ(mat1(1, 1), 6.0);
    EXPECT_EQ(mat1(1, 2), 7.0);
    EXPECT_EQ(mat1(2, 2), 9.0);
    EXPECT_EQ(mat1(0, 0), 1.0);
    mat1.print();
    print(mat1.row(1)); // Should print: [5,6,7]
    print(mat1.col(2)); // Should print: [7,9]
    EXPECT_TRUE((mat1.row(1)==il{5,6,7}));
    EXPECT_TRUE((mat1.col(2)==il{7,9}));
}
//,std::ranges::viewable_range V
template <std::ranges::range Range1,std::ranges::range Range2> auto Dot(const Range1& a, const Range2& b)
{
    double dot = 0.0;
    assert(a.size() == b.size() && "Ranges must be of the same size for dot product");
    for (auto [ia,ib]:std::views::zip(a,b)) {
        dot += ia*ib;
    }
    return dot;
}
template <class T, std::ranges::range Range> T Dot(const VectorView<Range>& a, const Vector<T>& b)
{
    assert(a.size() <= b.size() && "VectorView size exceeds Vector size");
    return Dot(a.range, b.view(a.indices).range);
}
template <class T, std::ranges::range Range> T Dot(const Vector<T>& a,const VectorView<Range>& b )
{
    assert(b.size() <= a.size() && "VectorView size exceeds Vector size");
    return Dot(a.view(b.indices).range, b.range);
}

template <typename T, class S> auto operator*(const Matrix<T,S>& m, const Vector<T>& v)
{
    auto indices = std::views::iota(size_t(0), v.size());
    auto dot=[m,v](size_t i) {return Dot(m.row(i),v);};
    auto mv=indices | std::views::transform(dot);
    return VectorView(std::move(mv), indices);
}

TEST_F(MatrixTests, DotProducts)
{
    TriDiagonalMatrix A{{1,2,0,0,0},
                        {5,6,7,0,0},
                        {0,8,9,10,0,0},
                        {0,0,11,12,13},
                        {0,0,0,14,15}};
    Vector<double> v1{1,2,3,4,5};
    EXPECT_EQ(Dot(A.row(0),v1), 1*1 + 2*2); 
    EXPECT_EQ(Dot(A.row(1),v1), 5*1 + 6*2 + 7*3); 
    EXPECT_EQ(Dot(A.row(2),v1), 8*2 + 9*3 + 10*4); 
    EXPECT_EQ(Dot(A.row(3),v1), 11*3 + 12*4 + 13*5); 
    EXPECT_EQ(Dot(A.row(4),v1), 14*4 + 15*5); 
    EXPECT_EQ(Dot(A.col(0),v1), 1*1 + 5*2); 
    EXPECT_EQ(Dot(A.col(1),v1), 2*1 + 6*2 + 8*3); 
    EXPECT_EQ(Dot(A.col(2),v1), 7*2 + 9*3 + 11*4); 
    EXPECT_EQ(Dot(A.col(3),v1), 10*3 + 12*4 + 14*5); 
    EXPECT_EQ(Dot(A.col(4),v1), 13*4 + 15*5); 
    
    EXPECT_EQ(Dot(v1,A.row(0)), 1*1 + 2*2); 
    EXPECT_EQ(Dot(v1,A.row(1)), 5*1 + 6*2 + 7*3); 
    EXPECT_EQ(Dot(v1,A.row(2)), 8*2 + 9*3 + 10*4); 
    EXPECT_EQ(Dot(v1,A.row(3)), 11*3 + 12*4 + 13*5); 
    EXPECT_EQ(Dot(v1,A.row(4)), 14*4 + 15*5); 
    EXPECT_EQ(Dot(v1,A.col(0)), 1*1 + 5*2); 
    EXPECT_EQ(Dot(v1,A.col(1)), 2*1 + 6*2 + 8*3); 
    EXPECT_EQ(Dot(v1,A.col(2)), 7*2 + 9*3 + 11*4); 
    EXPECT_EQ(Dot(v1,A.col(3)), 10*3 + 12*4 + 14*5); 
    EXPECT_EQ(Dot(v1,A.col(4)), 13*4 + 15*5);
 
    {
        Vector<double> Av(A * v1); // Matrix-vector multiplication
        EXPECT_EQ(Av(0), 1*1 + 2*2);
        EXPECT_EQ(Av(1), 5*1 + 6*2 + 7*3);
        EXPECT_EQ(Av(2), 8*2 + 9*3 + 10*4);
        EXPECT_EQ(Av(3), 11*3 + 12*4 + 13*5);
        EXPECT_EQ(Av(4), 14*4 + 15*5);
    }
    {
        Vector<double> Av=A * v1; // Matrix-vector multiplication
        EXPECT_EQ(Av(0), 1*1 + 2*2);
        EXPECT_EQ(Av(1), 5*1 + 6*2 + 7*3);
        EXPECT_EQ(Av(2), 8*2 + 9*3 + 10*4);
        EXPECT_EQ(Av(3), 11*3 + 12*4 + 13*5);
        EXPECT_EQ(Av(4), 14*4 + 15*5);
    }
}