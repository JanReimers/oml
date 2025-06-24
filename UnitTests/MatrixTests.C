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

struct intersection
{
    typedef std::ranges::iota_view<size_t,size_t> iota_view;
    intersection(const iota_view& a, const iota_view& b)
    {
        size_t i0=std::max(a.front(), b.front());
        size_t i1=std::min(a.back (), b.back ())+1;
        if (i0>i1) i1=i0;
        indices=std::ranges::iota_view(i0,i1); //new intersection range
        drop1 = a.front() < i0 ? i0 - a.front() : 0; // how many to drop from a
        drop2 = b.front() < i0 ? i0 - b.front() : 0; // how many to drop from b

    }

    iota_view indices;
    size_t drop1,drop2;
};

template <std::ranges::range Range1, std::ranges::range Range2> auto operator*(const VectorView<Range1>& a, const VectorView<Range2>& b)
{
    intersection inter(a.indices,b.indices);
    auto va=a.range | std::views::drop(inter.drop1) | std::views::take(inter.indices.size());
    auto vb=b.range | std::views::drop(inter.drop2) | std::views::take(inter.indices.size());
    return inner_product(va,vb);
}


template <class T, std::ranges::range Range> T operator*(const VectorView<Range>& a, const Vector<T>& b)
{
    assert(a.size() <= b.size() && "VectorView size exceeds Vector size");
    return a*b.view();
    // return inner_product(a.range, b.view(a.indices).range);
}
template <class T, std::ranges::range Range> T operator*(const Vector<T>& a,const VectorView<Range>& b )
{
    assert(b.size() <= a.size() && "VectorView size exceeds Vector size");
    return a.view()*b;
    // return inner_product(a.view(b.indices).range, b.range);
}
template <class T> T operator*(const Vector<T>& a,const Vector<T>& b )
{
    assert(b.size() == a.size() && "Vectors a and b size mismatch for dot product");
    return a.view()*b.view();
    // return inner_product(a.view(b.indices).range, b.range);
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
//     // assert(a.size() != 0 && b.size() != 0 && "Matrices must not be empty for multiplication");
//     assert(a.subscriptor.nc == b.subscriptor.nr && "Matrix dimensions do not match for multiplication");
//     auto indices = std::views::iota(size_t(0), a.subscriptor.nr);
//     auto ab=indices | std::views::transform([a,b](size_t i) {return a.row(i) * b;}); //uses op*(const Vector<T>& v,const Matrix<T,S>& m)
//     return Matrix<T,S>(std::move(ab), a.subscriptor.nr, b.subscriptor.nc);
// }


template <std::ranges::viewable_range R, std::ranges::viewable_range C, class S> class MatrixProductView
{
public:
    typedef std::ranges::range_value_t<R> T;
    MatrixProductView(const R& _rows, const C& _cols,S _subsciptor )
    : a_rows(_rows), b_cols(_cols), subscriptor(_subsciptor)
    {
        assert(nr()==subscriptor.nr);
        assert(nc()==subscriptor.nc);
    }
  
    size_t size() const { return  nr()*nc(); }
    size_t nr  () const { return std::ranges::size(a_rows); }
    size_t nc  () const { return std::ranges::size(b_cols); }

    auto operator()(size_t i, size_t j) const
    {
        // assert(subsciptor.is_stored(i,j) && "Index out of range for MatrixView");
        return a_rows[i]*b_cols[j]; //VectorView*VectorView
    }

    auto rows() const
    {
        auto outerp=std::views::cartesian_product(a_rows,b_cols) | std::views::transform([](auto tuple) {return get<0>(tuple) * get<1>(tuple);}); // uses op*(const Vector<T>& v,const Matrix<T,S>& m)
        return outerp | std::views::chunk(nc()) | std::views::transform([](auto chunk) {return VectorView(chunk);});
    }

    auto cols() const
    {
        auto outerp=std::views::cartesian_product(rows,cols) | std::views::transform([](auto tuple) {return get<0>(tuple) * get<1>(tuple);}); // uses op*(const Vector<T>& v,const Matrix<T,S>& m)
        return  std::views::iota(size_t(0), nc()) | std::views::transform
            ([outerp,this](size_t j) {return outerp | std::views::drop(j) | std::views::stride(nc());});
    }
    
// private:
    R a_rows; //a as a range fo rows.
    C b_cols; //b as a range of cols.
    S subscriptor; //packing for the product.
};


template <typename T, class S> auto operator*(const Matrix<T,S>& a,const Matrix<T,S>& b)
{
    assert(a.subscriptor.nc == b.subscriptor.nr && "Matrix dimensions do not match for multiplication");
    return MatrixProductView(a.rows(),b.cols(),a.subscriptor);
}

template <std::ranges::viewable_range R, std::ranges::viewable_range C,typename T, class S> 
auto operator*(const MatrixProductView<R,C,S>& a,const Matrix<T,S>& b)
{
    assert(a.nc() == b.subscriptor.nr && "Matrix dimensions do not match for multiplication");
    return MatrixProductView(a.rows(),b.cols(),a.subscriptor);
}

// template <std::ranges::viewable_range R1, std::ranges::viewable_range C1, std::ranges::viewable_range R2, std::ranges::viewable_range C2> 
// auto operator*(const MatrixView<R1,C1>& a, const MatrixView<R2,C2>& b)
// {
//     assert(a.nc()==b.nr() && "Matrix dimension mismatch for op*(M,M)");
//     auto outerp=std::views::cartesian_product(a.rows,b.cols) | std::views::transform([](auto tuple) {return get<0>(tuple) * get<1>(tuple);}); // uses op*(const Vector<T>& v,const Matrix<T,S>& m)
//     auto rows = outerp | std::views::chunk(b.nc()) | std::views::transform([](auto chunk) {return VectorView(chunk);});
//     auto cols = std::views::iota(size_t(0), b.nc()) | std::views::transform
//         ([outerp,b](size_t j) {return outerp | std::views::drop(j) | std::views::stride(b.nc());}
//     );
//     return MatrixView(rows, cols); // returns a MatrixView with the resulting rows and columns
// }



TEST_F(MatrixTests, DotProducts)
{
    Vector<double> v1{1,2,3,4,5};
    EXPECT_EQ(v1*v1,1*1+2*2+3*3+4*4+5*5);
    TriDiagonalMatrix A{{1,2,0,0,0},
                        {5,6,7,0,0},
                        {0,8,9,10,0,0},
                        {0,0,11,12,13},
                        {0,0,0,14,15}};
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
  
    intersection inter(i1,i2);
    auto v1_intersection=v1.range | std::views::drop(inter.drop1) | std::views::take(inter.indices.size());
    auto v2_intersection=v2.range | std::views::drop(inter.drop2) | std::views::take(inter.indices.size());
    int dot=inner_product(v1_intersection,v2_intersection);
    EXPECT_EQ(dot, 4*4 + 5*5 + 6*6 + 7*7 + 8*8); 
}

TEST_F(MatrixTests, ViewDotView2)
{
    std::valarray<int> v1(15),v2(13);
    std::ranges::iota(v1,0); // 0 1 2 3 4 5 6 7 8 9 10 11 12 13 14
    std::ranges::iota(v2,4); // 4 5 6 7 8 9 10 11 12 13 14 15 16
    auto i1=std::views::iota(3ul, 8ul+1ul); // these constant must be size_t = ul and not int which is the default.
    auto i2=std::views::iota(size_t(1), size_t(10+1)); //another way to do it.
    auto vv1=VectorView(v1 | std::views::drop(i1.front()) | std::views::take(i1.size()),i1); //     3 4 5 6 7 8 
    auto vv2=VectorView(v2 | std::views::drop(i2.front()) | std::views::take(i2.size()),i2); // 5 6 7 8 9 10 11 12 13 14
    EXPECT_EQ(vv1*vv2,3*7+4*8+5*9+6*10+7*11+8*12);
}

// TEST_F(MatrixTests, MatriView)
// {
//     TriDiagonalMatrix A{{1,2,0,0,0},
//                         {5,6,7,0,0},
//                         {0,8,9,10,0,0},
//                         {0,0,11,12,13},
//                         {0,0,0,14,15}};
//     auto rows=A.rows();
//     auto cols=A.cols();
//     auto view=A.view();
//     EXPECT_EQ(view.nr(), A.subscriptor.nr);
//     EXPECT_EQ(view.nc(), A.subscriptor.nc);
//     print2D(rows);
//     print2D(cols);
//     std::cout << "------------------------------------------------" << std::endl;
//     auto AA=view*view;
//     print2D(AA.rows);
//     print2D(AA.cols);
// }

TEST_F(MatrixTests, MatriProductView)
{
     TriDiagonalMatrix A{{1,2,0,0,0},
                        {5,6,7,0,0},
                        {0,8,9,10,0,0},
                        {0,0,11,12,13},
                        {0,0,0,14,15}};
    auto mpv=A*A; //MatriProductView
    EXPECT_EQ(mpv(0,0),11);
    EXPECT_EQ(mpv(0,1),14);
    EXPECT_EQ(mpv(0,2),14);
    EXPECT_EQ(mpv(1,0),35);
    EXPECT_EQ(mpv(1,1),102);
    EXPECT_EQ(mpv(1,2),105);
    EXPECT_EQ(mpv(1,3),70);
    EXPECT_EQ(mpv(2,0),40);
    EXPECT_EQ(mpv(2,1),120);
    EXPECT_EQ(mpv(2,2),247);
    EXPECT_EQ(mpv(2,3),210);
    EXPECT_EQ(mpv(2,4),130);
    EXPECT_EQ(mpv(3,1),88);
    EXPECT_EQ(mpv(3,2),231);
    EXPECT_EQ(mpv(3,3),436);
    EXPECT_EQ(mpv(3,4),351);
    EXPECT_EQ(mpv(4,2),154);
    EXPECT_EQ(mpv(4,3),378);
    EXPECT_EQ(mpv(4,4),407);

    auto mpvAAA=A*A*A;
    EXPECT_EQ(mpvAAA(0,0),81);
}