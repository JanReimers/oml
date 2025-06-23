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
template <std::ranges::range Range1,std::ranges::range Range2> auto inner_product(const Range1& a, const Range2& b)
{
    assert(a.size() == b.size() && "Ranges must be of the same size for dot product");
    std::ranges::range_value_t<Range1> dot(0); //
    for (auto [ia,ib]:std::views::zip(a,b)) dot += ia*ib;
    return dot;
} 

template <class T, std::ranges::range Range> T operator*(const VectorView<Range>& a, const Vector<T>& b)
{
    assert(a.size() <= b.size() && "VectorView size exceeds Vector size");
    return inner_product(a.range, b.view(a.indices).range);
}
template <class T, std::ranges::range Range> T operator*(const Vector<T>& a,const VectorView<Range>& b )
{
    assert(b.size() <= a.size() && "VectorView size exceeds Vector size");
    return inner_product(a.view(b.indices).range, b.range);
}

template <typename T, class S> auto operator*(const Matrix<T,S>& m, const Vector<T>& v)
{
    auto indices = std::views::iota(size_t(0), v.size());
    auto mv=indices | std::views::transform([m,v](size_t i) {return m.row(i) * v;});
    return VectorView(std::move(mv), indices);
}
template <typename T, class S> auto operator*(const Vector<T>& v,const Matrix<T,S>& m)
{
    auto indices = std::views::iota(size_t(0), v.size());
    auto vm=indices | std::views::transform([m,v](size_t j) {return v * m.col(j);});
    return VectorView(std::move(vm), indices);
}
template <typename T, class S,std::ranges::range R> auto operator*(const VectorView<R>& v,const Matrix<T,S>& m)
{
    auto indices = std::views::iota(size_t(0), v.size());
    auto vm=indices | std::views::transform([m,v](size_t j) {return v * m.col(j);});
    return VectorView(std::move(vm), indices);
}
// template <typename T, class S> auto operator*(const Matrix<T,S>& a,const Matrix<T,S>& b)
// {
//     assert(a.size() != 0 && b.size() != 0 && "Matrices must not be empty for multiplication");
//     assert(a.subscriptor.nc == b.subscriptor.nr && "Matrix dimensions do not match for multiplication");
//     auto indices = std::views::iota(size_t(0), a.subscriptor.nr);
//     auto ab=indices | std::views::transform([a,b](size_t i) {return a.row(i) * b;}); //uses op*(const Vector<T>& v,const Matrix<T,S>& m)
//     return Matrix<T,S>(std::move(ab), a.subscriptor.nr, b.subscriptor.nc);
// }


TEST_F(MatrixTests, DotProducts)
{
    TriDiagonalMatrix A{{1,2,0,0,0},
                        {5,6,7,0,0},
                        {0,8,9,10,0,0},
                        {0,0,11,12,13},
                        {0,0,0,14,15}};
    Vector<double> v1{1,2,3,4,5};
    EXPECT_EQ(A.row(0)*v1, 1*1 + 2*2); 
    EXPECT_EQ(A.row(1)*v1, 5*1 + 6*2 + 7*3); 
    EXPECT_EQ(A.row(2)*v1, 8*2 + 9*3 + 10*4); 
    EXPECT_EQ(A.row(3)*v1, 11*3 + 12*4 + 13*5); 
    EXPECT_EQ(A.row(4)*v1, 14*4 + 15*5); 
    EXPECT_EQ(A.col(0)*v1, 1*1 + 5*2); 
    EXPECT_EQ(A.col(1)*v1, 2*1 + 6*2 + 8*3); 
    EXPECT_EQ(A.col(2)*v1, 7*2 + 9*3 + 11*4); 
    EXPECT_EQ(A.col(3)*v1, 10*3 + 12*4 + 14*5); 
    EXPECT_EQ(A.col(4)*v1, 13*4 + 15*5); 
    
    EXPECT_EQ(v1*A.row(0), 1*1 + 2*2); 
    EXPECT_EQ(v1*A.row(1), 5*1 + 6*2 + 7*3); 
    EXPECT_EQ(v1*A.row(2), 8*2 + 9*3 + 10*4); 
    EXPECT_EQ(v1*A.row(3), 11*3 + 12*4 + 13*5); 
    EXPECT_EQ(v1*A.row(4), 14*4 + 15*5); 
    EXPECT_EQ(v1*A.col(0), 1*1 + 5*2); 
    EXPECT_EQ(v1*A.col(1), 2*1 + 6*2 + 8*3); 
    EXPECT_EQ(v1*A.col(2), 7*2 + 9*3 + 11*4); 
    EXPECT_EQ(v1*A.col(3), 10*3 + 12*4 + 14*5); 
    EXPECT_EQ(v1*A.col(4), 13*4 + 15*5);
 
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
    {
        Vector<double> vA = v1 * A; // Vector-Matrix multiplication
        EXPECT_EQ(vA(0), 1*1 + 5*2 + 0*3 + 0*4 + 0*5);
        EXPECT_EQ(vA(1), 2*1 + 6*2 + 8*3 + 0*4 + 0*5);
        EXPECT_EQ(vA(2), 0*1 + 7*2 + 9*3 + 11*4 + 0*5);
        EXPECT_EQ(vA(3), 0*1 + 0*2 + 10*3 + 12*4 + 14*5);
        EXPECT_EQ(vA(4), 0*1 + 0*2 + 0*3 + 13*4 + 15*5);
        print(vA); // Should print: [11,38,85,148,127]
    }
    EXPECT_EQ(v1*A*v1,11+2*38+3*85+4*148+5*127);
    
    // v1.view()*A; need to figure out how to handle indices for inner product with views.
}

#include <ranges>
#include <vector>
// #include<algorithm>
#include <numeric>
TEST_F(MatrixTests, ViewDotView)
{
    std::valarray<int> v(15);
    std::ranges::iota(v,0); // need to include <numeric> for iota!
    auto i1=std::views::iota(size_t(3), size_t(8+1));
    auto i2=std::views::iota(size_t(4), size_t(10+1));
    auto v1v=v | std::views::drop(i1.front()) | std::views::take(i1.size());
    auto v2v=v | std::views::drop(i2.front()) | std::views::take(i2.size());
    auto v1=VectorView(v1v,i1); // Full view of the valarray
    auto v2=VectorView(v2v,i2); // Full view of the valarray
    print(v);
    print(v1v);
    print(v2v);
    print(v1);
    print(v2);
    typedef std::ranges::iota_view<size_t,size_t> iota_view;
    struct intersection
    {
        intersection(const iota_view& a, const iota_view& b)
        {
            size_t i0=std::max(a.front(), b.front());
            size_t i1=std::min(a.back (), b.back ());
            indices=std::ranges::iota_view(i0,i1+1); //new intersection range
            drop1 = a.front() < i0 ? i0 - a.front() : 0; // how many to drop from a
            drop2 = b.front() < i0 ? i0 - b.front() : 0; // how many to drop from b

        }
        iota_view indices;
        size_t drop1,drop2;
    };
    intersection inter(i1,i2);
    auto v1_intersection=v1.range | std::views::drop(inter.drop1) | std::views::take(inter.indices.size());
    auto v2_intersection=v2.range | std::views::drop(inter.drop2) | std::views::take(inter.indices.size());
    print(v1_intersection);
    print(v2_intersection);
    int dot=inner_product(v1_intersection,v2_intersection);
    EXPECT_EQ(dot, 4*4 + 5*5 + 6*6 + 7*7 + 8*8); 
}