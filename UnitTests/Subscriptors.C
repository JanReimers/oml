// File: Subscriptors.C  Experiment with various matrix packings.

// Notes:
//   1) how to load an upper triangular from an initializer list.
//
#include "gtest/gtest.h"
class SubscriptorsTests : public ::testing::Test
{
public:
    SubscriptorsTests() = default;
    ~SubscriptorsTests() override = default;
};

#include <ranges>

class FullSubsciptor
{
public:
    FullSubsciptor(size_t nrows, size_t ncols)
        : nr(nrows), nc(ncols) {};
    bool is_stored(size_t i, size_t j) const
    {
        return i>=0 && i<nr && j>=0 && j<nc; // Full matrix, all elements are stored
    }
    size_t offset(size_t i, size_t j) const
    {
        return i * nc + j;
    }
    size_t size() const
    {
        return nr * nc; // Total number of elements
    }
    size_t stored_row_size(size_t) const //Don't need the row index.
    {
        return nc; // Each row has ncols elements
    }
    template <std::ranges::range R> auto view(R& r) const
    {
        return std::views::chunk(r,stored_row_size(0));
    }
    size_t nr,nc;
};

class UpperTriangularSubsciptor
{
public:
    UpperTriangularSubsciptor(size_t nrows, size_t ncols)
        : nr(nrows), nc(ncols) {};
    bool is_stored(size_t i, size_t j) const
    {
        return i <= j; 
    }
    size_t offset(size_t i, size_t j) const
    {
        assert(is_stored(i,j));
        return i * (2*nc-i-1) / 2 + j; // Upper triangular matrix        
    }
    size_t size() const
    {
        return nc<=nr ? nc*(nc+1)/2 : nr*(nr+1)/2 +(nc-nr)*nr; // Total number of elements
    }
    size_t stored_row_size(size_t row_index) const
    {
        assert(row_index < nr);
        if (row_index >= nc) return 0; // No elements stored in this row
        return nc-row_index + ((nc>nr) ? nc-nr : 0); // Each row has ncols elements
    }
    template <std::ranges::range R> auto view(R& r) const
    {
        auto row = [r,this](int i) mutable 
        {
            return r | std::views::drop(offset(i,i)) | std::views::take(stored_row_size(i));
        };
        return std::views::iota(0,(int)nr) | std::views::transform(row);
    }
    size_t nr,nc;
};

class DiagonalSubsciptor
{
    public:
     DiagonalSubsciptor(size_t nrows, size_t ncols)
        : nr(nrows), nc(ncols) {};
    bool is_stored(size_t i, size_t j) const
    {
        return i == j; 
    }
    size_t offset(size_t i, size_t j) const
    {
        assert(is_stored(i,j));
        return i; // Diagonal matrix
    }
    size_t size() const
    {
        return std::min(nr,nc); // Total number of elements
    }
    size_t stored_row_size(size_t) const
    {
        return 1;
    }
    template <std::ranges::range R> auto view(R& r) const
    {
        auto row = [r,this](int i) mutable 
        {
            return r | std::views::drop(offset(i,i)) | std::views::take(stored_row_size(i));
        };
        return std::views::iota(0,(int)nr) | std::views::transform(row);
    }
    size_t nr,nc;
};

class TriDiagonalSubsciptor
{
    public:
     TriDiagonalSubsciptor(size_t nrows, size_t ncols)
        : nr(nrows), nc(ncols) {};
    bool is_stored(size_t i, size_t j) const
    {
        return i<=j+1 && j<=i+1; // Main diagonal and two adjacent diagonals
    }
    size_t offset(size_t i, size_t j) const
    {
        assert(is_stored(i,j));
        return 2*i+j; // Tri-diagonal matrix
    }
    size_t size() const
    {
        return nc<=nr ? nc*(nc+1)/2 : nr*(nr+1)/2 +(nc-nr)*nr; // Total number of elements
    }
    size_t stored_row_size(size_t row_index) const
    {
        return (row_index==0 || row_index==nr-1) ? 2 : 3;
    }
    template <std::ranges::range R> auto view(R& r) const
    {
        auto row = [r,this](int i) mutable 
        {
            size_t j=i>0 ? i-1 : 0;
            return r | std::views::drop(offset(i,j)) | std::views::take(stored_row_size(i));
        };
        return std::views::iota(0,(int)nr) | std::views::transform(row);
    }
    size_t nr,nc;
};

class BandedSubsciptor
{
    public:
    BandedSubsciptor(size_t nrows, size_t ncols, size_t _k)
        :  nr(nrows), nc(ncols), k(_k) {}

    bool is_stored(size_t i, size_t j) const
    {
        return (i<=j+k && j<=i+k);
    }

    size_t offset(size_t i, size_t j) const
    {
        assert(is_stored(i,j));
        size_t dj=j;
        if (i>k) dj+=i;
        return (i * ( k + 1)) + dj; // Banded matrix
    }
    size_t size() const
    {
        return (nr * (2*k + 1)) - k * (k + 1) / 2; // Total number of elements
    
    }
    size_t stored_row_size(size_t row_index) const
    {
        if (row_index < k) return row_index + 1; // First k rows have increasing size
        if (row_index >= nr - k) return nc - row_index + k; // Last k rows have decreasing size
        return 2 * k + 1; // Middle rows have full band width
    }
    template <std::ranges::range R> auto view(R& r) const
    {
        auto row = [r,this](int i) mutable 
        {
            return r | std::views::drop(offset(i,i)) | std::views::take(stored_row_size(i));
        };
        return std::views::iota(0,(int)nr) | std::views::transform(row);
    }
    
    size_t nr,nc,k;

};

#include <valarray>
#include <cassert>
template <typename T, class S> class Matrix
{
    public:
    Matrix(size_t nr, size_t nc) : subscriptor(nc,nr), data(subscriptor.size()) {};
    Matrix(std::initializer_list<std::initializer_list<T>> init)
        : subscriptor(init.size(), init.begin()->size()), data(subscriptor.size())
    {
        load(init);
    }
    void load(std::initializer_list<std::initializer_list<T>> init)
    {
        assert(init.size() == subscriptor.nr && "Initializer list size does not match subscriptor row count");
        assert(init.begin()->size() == subscriptor.nc && "Initializer list row size does not match subscriptor column count");
        size_t i = 0;
        for (const auto& row : init)
        {
            size_t j = 0;
            for (const auto& val : row)
            {
                if (subscriptor.is_stored(i, j))
                    data[subscriptor.offset(i, j)] = val;
                else
                    assert(val==0.0);
                ++j;
            }
            ++i;
        }
    }
    Matrix(size_t nr, size_t nc, size_t k) : subscriptor(nc,nr,k), data(subscriptor.size()) {};
    Matrix(std::initializer_list<std::initializer_list<T>> init, size_t k)
    : subscriptor(init.size(), init.begin()->size(),k), data(subscriptor.size())
    {
        load(init);
    }
    T operator()(size_t i, size_t j) const
    {
        assert(subscriptor.is_stored(i,j));
        return data[subscriptor.offset(i,j)];
    }
    T& operator()(size_t i, size_t j)
    {
        assert(subscriptor.is_stored(i,j));
        return data[subscriptor.offset(i,j)];
    }
    size_t size() const
    {
        return subscriptor.size();
    }
    auto begin()       { return std::begin(data); }
    auto end  ()       { return std::end  (data); }
    auto begin() const { return std::begin(data); }
    auto end  () const { return std::end  (data); }
    auto view () 
    {
        return subscriptor.view(data);
    }

    void print() const
    {
        for (size_t i = 0; i < subscriptor.nr; ++i)
        {
            for (size_t j = 0; j < subscriptor.nc; ++j)
            {
                if (subscriptor.is_stored(i, j))
                    std::cout << (*this)(i, j) << " ";
                else
                    std::cout << "0 "; // Print 0 for non-stored elements
            }
            std::cout << std::endl;
        }   
    }
private:
    S subscriptor;
    std::valarray<T> data;
};

template <std::ranges::range Range> void print(Range v)
{
    for (auto element : v) {
        std::cout << element << " ";
    }
    std::cout << "\n";
}
auto print2D_3=[](auto rng){for(auto r:rng)print(r);};

TEST_F(SubscriptorsTests, FullSubscriptor)
{
    FullSubsciptor sub(3, 4);
    EXPECT_TRUE(sub.is_stored(0, 0));
    EXPECT_TRUE(sub.is_stored(2, 3));
    EXPECT_FALSE(sub.is_stored(3, 0));
    EXPECT_FALSE(sub.is_stored(0, 4));
    EXPECT_EQ(sub.offset(1, 2), 6);
    EXPECT_EQ(sub.size(), 12);
    Matrix<double, FullSubsciptor> mat(3, 4);
    mat(1, 2) = 5.0;
    EXPECT_EQ(mat(1, 2), 5.0);
    Matrix<double, FullSubsciptor> mat1{{1,2,3},
                                       {4,5,6},
                                       {7,8,9},
                                       {10,11,12}};
    EXPECT_EQ(mat1(1, 2), 6.0);
    EXPECT_EQ(mat1(2, 1), 8.0);
    EXPECT_EQ(mat1(0, 0), 1.0);
    EXPECT_EQ(mat1(0, 1), 2.0);
    EXPECT_EQ(mat1(0, 2), 3.0);
    mat1.print();
    print2D_3(mat1.view());
}
TEST_F(SubscriptorsTests, UpperTriangularSubscriptor)
{
    {
    UpperTriangularSubsciptor sub(3,3);
    EXPECT_TRUE(sub.is_stored(0, 0));
    EXPECT_TRUE(sub.is_stored(1, 1));
    EXPECT_TRUE(sub.is_stored(2, 2));
    EXPECT_TRUE(sub.is_stored(1, 2));
    EXPECT_FALSE(sub.is_stored(2, 1));
    }
    {
    Matrix<double, UpperTriangularSubsciptor> A(3, 4);
    A(1, 2) = 5.0;
    EXPECT_EQ(A(1, 2), 5.0);
    }
    {
        Matrix<double, UpperTriangularSubsciptor> A{{1,2,3,4},
                                       {0,5,6,7},
                                       {0,0,8,9},
                                       {0,0,0,10}};
        A.print();
        print2D_3(A.view());
    }

}
TEST_F(SubscriptorsTests, DiagonalSubscriptor)
{
    {
        DiagonalSubsciptor sub(3,4);
        EXPECT_TRUE(sub.is_stored(0, 0));
        EXPECT_TRUE(sub.is_stored(1, 1));
        EXPECT_TRUE(sub.is_stored(2, 2));
        EXPECT_FALSE(sub.is_stored(1, 2));
        EXPECT_EQ(sub.offset(1, 1), 1);
    }
    
    {
        Matrix<double, DiagonalSubsciptor> mat(3, 4);
        mat(1, 1) = 5.0;
        EXPECT_EQ(mat(1, 1), 5.0);
    }
    {
        Matrix<double, DiagonalSubsciptor> A{{1,0,0,0},
                                       {0,5,0,0},
                                       {0,0,8,0},
                                       {0,0,0,10}};
        A.print();
        print2D_3(A.view());
    }
}
TEST_F(SubscriptorsTests, TriDiagonalSubscriptor)
{
    {
        TriDiagonalSubsciptor sub(3,4);
        EXPECT_TRUE(sub.is_stored(0, 0));
        EXPECT_TRUE(sub.is_stored(1, 1));
        EXPECT_TRUE(sub.is_stored(2, 2));
        EXPECT_TRUE(sub.is_stored(1, 0));
        EXPECT_TRUE(sub.is_stored(1, 2));
        EXPECT_FALSE(sub.is_stored(2, 0));
        EXPECT_EQ(sub.offset(1, 2), 4);
    }
    {
        Matrix<double, TriDiagonalSubsciptor> mat(3, 4);
        mat(1, 2) = 5.0;
        EXPECT_EQ(mat(1, 2), 5.0);
    }
    {
        Matrix<double, TriDiagonalSubsciptor> A{{1,2,0,0},
                                       {3,4,5,0},
                                       {0,6,7,8},
                                       {0,0,9,10}};
        A.print();
        print2D_3(A.view());
    }
}
// TEST_F(SubscriptorsTests, BandedSubscriptor)
// {
//     {
//         BandedSubsciptor sub(5, 5, 2);
//         EXPECT_TRUE(sub.is_stored(0, 0));
//         EXPECT_TRUE(sub.is_stored(1, 0));
//         EXPECT_TRUE(sub.is_stored(1, 1));
//         EXPECT_TRUE(sub.is_stored(2, 1));
//         EXPECT_TRUE(sub.is_stored(2, 2));
//         EXPECT_TRUE(sub.is_stored(3, 2));
//         EXPECT_TRUE(sub.is_stored(3, 3));
//         EXPECT_TRUE(sub.is_stored(4, 3));
//         EXPECT_FALSE(sub.is_stored(0, 4));
//         EXPECT_FALSE(sub.is_stored(1, 4));
//         EXPECT_FALSE(sub.is_stored(2, 5));
//
//         EXPECT_EQ(sub.offset(0, 0), 0);
//         EXPECT_EQ(sub.offset(0, 1), 1);
//         EXPECT_EQ(sub.offset(0, 2), 2);
//         EXPECT_EQ(sub.offset(1, 0), 3);
//         EXPECT_EQ(sub.offset(1, 1), 4);
//         EXPECT_EQ(sub.offset(1, 2), 5);
//         EXPECT_EQ(sub.offset(1, 3), 6);
//         EXPECT_EQ(sub.offset(2, 0), 7);
//         EXPECT_EQ(sub.offset(2, 1), 8);
//         EXPECT_EQ(sub.offset(2, 2), 9);
//         EXPECT_EQ(sub.offset(2, 3), 10);
//         EXPECT_EQ(sub.offset(2, 4), 11);
//     }
//     {
//         Matrix<double, BandedSubsciptor> mat(5, 5,2);
//         mat(1, 2) = 5.0;
//         EXPECT_EQ(mat(1, 2), 5.0);
//         mat(1, 3) = 5.1;
//         EXPECT_EQ(mat(1, 3), 5.1);
//     }
//     {
//         Matrix<double, BandedSubsciptor> A({{1,2,0,0,0},
//                                        {3,4,5,0,0},
//                                        {0,6,7,8,0},
//                                        {0,0,9,10,11},
//                                        {0,0,12,13,14}},2);
//         A.print();
//         print2D_3(A.view());
//     }
// }