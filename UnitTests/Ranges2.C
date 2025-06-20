// File: Ranges2.C  Experiment with C++-20 ranges.

#include "gtest/gtest.h"
#include <ranges>
#include <iostream>
#include <vector>
#include <algorithm>
#include <generator>
#include <valarray>
#include <cassert>

using std::cout;
using std::endl;

class RangesTests2 : public ::testing::Test
{
public:
};

template <std::ranges::range Range> void print(Range v)
{
    for (auto element : v) {
        std::cout << element << " ";
    }
    std::cout << "\n";
}
auto print2D_2=[](auto rng){for(auto r:rng)print(r);};

template <std::ranges::viewable_range R> class VectorView
{
public:
    VectorView(R&& r, size_t _n)
        : range(std::forward<R>(r)), n(_n)
    {}

    auto begin()       { return range.begin(); }
    auto end  ()       { return range.end  (); }
    auto begin() const { return range.begin(); }
    auto end  () const { return range.end  (); }
    // many ranges don't support random access with op[].
    size_t size() const { return  n; }
    R range; // The underlying range
    size_t n;
};
template <std::ranges::viewable_range R> class MatrixView
{
public:
    MatrixView(R&& r, size_t _n, size_t _nr, size_t _nc)
        : range(std::forward<R>(r)), n(_n), nr(_nr), nc(_nc)
    {}

    auto begin()       { return range.begin(); }
    auto end  ()       { return range.end  (); }
    auto begin() const { return range.begin(); }
    auto end  () const { return range.end  (); }
    // many ranges don't support random access with op[].
    size_t size() const { return  n; }
    R range; // The underlying range
    size_t n,nr,nc;
};

class Vector
{
public:
    Vector() : data(0) {};
   
    Vector(size_t n)
        : data(n)
    {}

    template <std::ranges::view V> Vector(VectorView<V>& view) //const won't work for views.
        : data(view.size())
    {
        assign_from(view);
    }
    template <std::ranges::view V> Vector(VectorView<V>&& view) //const won't work for views.
        : data(view.size())
    {
        assign_from(view);
    }
    template <std::ranges::view V> Vector& operator=(VectorView<V>&& view) //const won't work for views.
    {   
        if (data.size() == 0)
            data.resize(view.size());
        else
            assert(size()==view.size());

        assign_from(view);
        return *this;
    }

    double operator()(size_t i) const
    {
        assert(i < size());
        return data[i];
    }
    double& operator()(size_t i)
    {
        assert(i < size());
        return data[i];
    }
    size_t size() const { return data.size(); }
    auto  begin()       { return std::begin(data); }
    auto  end  ()       { return std::end  (data); }
    auto  begin() const { return std::begin(data); }
    auto  end  () const { return std::end  (data); }
private:
    template <std::ranges::view V> void assign_from(VectorView<V>& view)
    {
        size_t i=0;
        for (auto& r:view.range) data[i++] = r; //This should be where the lazy evaluation of all the chained views happens.
    }

    std::valarray<double> data;
};

class Matrix
{
public:
    Matrix(size_t nrows, size_t ncols)
        : nr(nrows), nc(ncols), data(nr*nc)
    {}
    Matrix(std::initializer_list<std::initializer_list<double>> init)
        : nr(init.size()), nc(init.begin()->size()), data(nr * nc)
    {
        size_t i = 0;
        for (const auto& row : init) {
            assert(row.size() == nc); // Ensure all rows have the same number of columns
            std::copy(row.begin(), row.end(), &data[i * nc]);
            ++i;
        }
    }
    template <std::ranges::view V> Matrix(MatrixView<V>&& view) //const won't work for views.
        : data(view.size())
    {
        assign_from(view);
    }

    size_t size() const { return data.size(); } // Total number of elements
    auto  begin()       { return std::begin(data); }
    auto  end  ()       { return std::end  (data); }
    auto  begin() const { return std::begin(data); }
    auto  end  () const { return std::end  (data); }
    double operator()(size_t i, size_t j) const
    {
        assert(i < nr && j < nc);
        return data[i * nc + j];
    }
    double& operator()(size_t i, size_t j)
    {
        assert(i < nr && j < nc);
        return data[i * nc + j];
    }

    auto row(size_t i)
    {
        assert(i < nr);
        auto v= data | std::views::chunk(nc) 
                     | std::views::drop(i) 
                     | std::views::take(1) 
                     | std::views::join; // Return a view of the ith row
        return VectorView<decltype(v)>(std::move(v),nc);
    }
    auto col(size_t j)
    {
        assert(j < nc);
        auto v= data | std::views::drop(j) 
                     | std::views::stride(nc) 
                     | std::views::take(nr); // Return a view of the jth column
        return VectorView<decltype(v)>(std::move(v),nr);
    }

    auto transpose1() 
    {
        namespace rs = std::ranges;
        namespace vs = std::ranges::views;
        auto get_column = [*this](int i) mutable 
        {
            return data | vs::drop(i) | vs::stride(nc);
        };
        auto m= vs::iota(0,(int)nc) | vs::transform(get_column);
        return MatrixView<decltype(m)>(std::move(m), size(), nc, nr);
    }

private:
    template <std::ranges::view V> void assign_from(MatrixView<V>& view)
    {
        size_t i=0;
        for (auto r:view.range | std::views::join) data[i++] = r; //This should be where the lazy evaluation of all the chained views happens.

    }
    size_t nr, nc;
    std::valarray<double> data;
};

TEST_F(RangesTests2, MatrixViewDemo)
{
    Matrix A{{1,2,3,4},{5,6,7,8},{9,10,11,12}};

    auto row_view = A.row(1); // Get the second row, returns a VectorView proxy
    print(row_view); // Should print: [5,6,7,8]
    EXPECT_EQ(row_view.size(), 4);

    auto col_view = A.col(2); // Get the third column returns a VectorView proxy
    print(col_view); // Should print: [3,7,11]
    EXPECT_EQ(col_view.size(), 3);
   

    Vector vrow1(row_view);
    print(vrow1);
    EXPECT_EQ(vrow1.size(), 4);
    EXPECT_EQ(vrow1(0), 5);
    EXPECT_EQ(vrow1(1), 6); 
    EXPECT_EQ(vrow1(2), 7);
    EXPECT_EQ(vrow1(3), 8);
    vrow1= A.row(1);
    print(vrow1);
    EXPECT_EQ(vrow1.size(), 4);
    EXPECT_EQ(vrow1(0), 5);
    EXPECT_EQ(vrow1(1), 6); 
    EXPECT_EQ(vrow1(2), 7);
    EXPECT_EQ(vrow1(3), 8);
    Vector vrow2;
    vrow2= A.row(1); //op= should resize
    EXPECT_EQ(vrow2.size(), 4);
    EXPECT_EQ(vrow2(0), 5);
    EXPECT_EQ(vrow2(1), 6); 
    EXPECT_EQ(vrow2(2), 7);
    EXPECT_EQ(vrow2(3), 8);

    Vector vcol=col_view; // Get the third column
    EXPECT_EQ(vcol.size(), 3);
    EXPECT_EQ(vcol(0), 3);
    EXPECT_EQ(vcol(1), 7);
    EXPECT_EQ(vcol(2), 11);
//
//  Now try and avoid the proxy intermetiadte
//
    Vector vr_direct(A.row(1)),vc_direct(A.col(2)); // Get row or column directly
    EXPECT_EQ(vr_direct.size(), 4);
    EXPECT_EQ(vr_direct(0), 5);
    EXPECT_EQ(vr_direct(1), 6); 
    EXPECT_EQ(vr_direct(2), 7);
    EXPECT_EQ(vr_direct(3), 8);
    EXPECT_EQ(vc_direct.size(), 3);
    EXPECT_EQ(vc_direct(0), 3);
    EXPECT_EQ(vc_direct(1), 7);
    EXPECT_EQ(vc_direct(2), 11);
}



TEST_F(RangesTests2, MatrixTransposeDemo)
{
    Matrix A{{1,2,3,4},{5,6,7,8},{9,10,11,12}};
    auto At=A.transpose1();
    print2D_2(At);
    Matrix At1(A.transpose1());
    print(At1);
}

