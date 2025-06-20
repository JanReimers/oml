// File: Ranges.C  Experiment with C++-20 ranges.
//  Following https://gieseanw.wordpress.com/2019/10/20/we-dont-need-no-stinking-expression-templates/

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

class RangesTests : public ::testing::Test
{
public:
};

TEST_F(RangesTests, BasicRange)
{
    // Example of using C++20 ranges
    std::vector<int> vec = {1, 2, 3, 4, 5};

    // Using ranges to filter and transform
    auto result = vec | std::views::filter([](int x) { return x % 2 == 0; })
                      | std::views::transform([](int x) { return x * x; });

    std::vector<int> squared_evens(result.begin(), result.end());
    // std::vector<int> squared_evens(std::from_range,result); c++-23 version.

    EXPECT_EQ(squared_evens.size(), 2);
    EXPECT_EQ(squared_evens[0], 4); // 2^2
    EXPECT_EQ(squared_evens[1], 16); // 4^2
}

// auto zip_with(auto op, auto&&... ranges)
// {
//     return std::views::zip(std::forward<decltype(ranges)>(ranges)...)
//            | std::views::transform([op](auto&&... args) {
//                  return (op(std::forward<decltype(args)>(args)...));
//              });
// }

template <std::ranges::range Range> void print(Range v)
{
    for (auto element : v) {
        std::cout << element << " ";
    }
    std::cout << "\n";
}

// Try out std::views::zip_transform
TEST_F(RangesTests, zip_transform_Demo1)
{
    std::vector<int> a = {1, 2, 3, 4, 5};
    std::vector<int> b = {10, 20, 30, 40, 50};
    auto plus = [](int x, int y) { return x + y; };
    auto p=std::views::zip_transform(plus,a,b);
    print(p);
    std::ranges::for_each(p, [](auto n) { std::cout << n << " "; });
    std::cout << std::endl;
    for (auto [pi,ai,bi] : std::views::zip(p,a,b)) {
        EXPECT_EQ(pi,ai+bi);
    }
  
}

// Try out std::views::zip_transform with singles.
TEST_F(RangesTests, zip_transform_Demo2)
{
    int left = 1;
    double right = 40.0;
    auto plus = [](int x, int y) { return x + y; };
    auto expr1=std::views::zip_transform(plus,std::ranges::single_view{left}, std::ranges::single_view{right});
    print(expr1);
    EXPECT_EQ(*expr1.begin(), 41);
    EXPECT_EQ(expr1[0], 41);
    short next_value = 1;
    auto expr2 = std::views::zip_transform(std::plus<>{}, expr1, std::ranges::single_view{next_value});
    print(expr2);
    EXPECT_EQ(*expr2.begin(), 42);
    EXPECT_EQ(expr2[0], 42);
}

TEST_F(RangesTests, ScaledVectorAddition)
{
    std::vector<int> left_vals{1, 2, 3};
    std::vector<int> right_vals{1, 1, 1};
    auto alpha = 4;
    auto scaled_rhs = right_vals | std::views::transform([&alpha](int val){return alpha * val;});
    auto scaled_addition = std::views::zip_transform(std::plus<>{}, left_vals, scaled_rhs);

    print(scaled_addition);
    for (auto [pi,ai,bi] : std::views::zip(scaled_addition,left_vals,right_vals)) {
        EXPECT_EQ(pi,ai+alpha*bi);
    }
    
}

// std::ranges::range is a concept (not a type).
template <std::ranges::range Range> double my_inner_product(const Range& a, const Range& b, double init)
{
    for (auto&& [ai, bi] : std::views::zip(a, b)) {
        init += ai * bi;
    }
    return init;
}

TEST_F(RangesTests, MV_Mutiply1)
{
    using Matrix = std::vector<std::vector<double>>;
    using Vector = std::vector<double>;

    Matrix A = {{1, 2, 3}, {4, 5, 6}};
    Vector x = {7,8, 9};
    auto product = A | std::views::transform([&x](const auto& row){
                            return my_inner_product(row, x, 0.0);
                        });

    EXPECT_EQ(product[0],50);
    EXPECT_EQ(product[1],122);
}

//Try with vallarray instead of vector.
TEST_F(RangesTests, MV_Mutiply2)
{
    using Matrix = std::valarray<std::valarray<double>>;
    using Vector = std::valarray<double>;

    Matrix A = {{1, 2, 3}, {4, 5, 6}};
    Vector x = {7,8, 9};
    auto product = A | std::views::transform([&x](const auto& row){
                            return my_inner_product(row, x, 0.0);
                        });

    EXPECT_EQ(product[0],50);
    EXPECT_EQ(product[1],122);
}


using Matrix = std::valarray<std::valarray<double>>;
using Vector = std::valarray<double>;
// Ax = B
auto operator*(const Matrix& A, const Vector& x)
{
    return A | std::views::transform
    ( [&x](const auto& row)
        {
            return my_inner_product(row, x, 0.0);
        }
    );
}

TEST_F(RangesTests, MV_Mutiply3)
{
    Matrix A = {{1, 2, 3}, {4, 5, 6}};
    Vector x = {7,8, 9};
    auto Ax=A*x;
    EXPECT_EQ(Ax[0],50);
    EXPECT_EQ(Ax[1],122);
    std::vector<double> Ax2(Ax.begin(), Ax.end());
    std::vector<double> Ax3(std::from_range,Ax); //c++-23
    std::vector<double> Ax4(std::from_range,A*x); //c++-23
}

// std::generator<double> as_generator(const Vector& vec)
// {
//     for(double val : vec)
//         co_yield val;
// }

std::vector<double> as_stdvector(std::generator<double> gen)
{
    std::vector<double> result;
    for (auto val : gen) {
        result.push_back(val);
    }
    return result;
}

std::generator<double> operator*(const Vector& vec, double factor)
{
    for(double val : vec)
        co_yield val * factor;
}
std::generator<double> operator*(double factor,const Vector& vec)
{
    for(double val : vec)
        co_yield factor * val;
}

std::generator<double> operator+( std::generator<double> lhs,  std::generator<double> rhs)
{
    auto ilhs=lhs.begin();
    for (auto irhs=rhs.begin();irhs!=rhs.end();irhs++,ilhs++ )
        co_yield *ilhs+*irhs;
}

TEST_F(RangesTests, GeneratorDemo)
{
    Vector v1 = {1, 2, 3};
    Vector v2 = {4, 5, 6};
    std::vector<double> result=as_stdvector(v1 * 2.0 + 3.0 * v2);
    
    print(result);
    EXPECT_EQ(result.size(), 3);
    EXPECT_EQ(result[0], 14); // 2*1 + 3*4
    EXPECT_EQ(result[1], 19); // 2*2 + 3*5
    EXPECT_EQ(result[2], 24); // 2*3 + 3*6
}

// Got stuck here.  We something like a 2D generator std::generator<std::generator<double>> does not work.
// std::ranges::view<std::generator<double>> transposed(const Matrix& mat)
// {
//     assert(mat.size()!=0);
//     size_t nrows=mat.size(), ncols = mat[0].size();
//     for(size_t j = 0; j < ncols; ++j)
//         for(size_t i = 0; i < nrows; ++i)
//             co_yield mat[i][j];
// }


// TEST_F(RangesTests, TransposeDemo)
// {
//     Matrix A = {{1, 2, 3}, {4, 5, 6}};
//     auto At = transposed(A);
//     EXPECT_EQ(At[0][1],5); // 2 rows * 3 cols
// }

// // x*A = B

// Now look at https://mmore500.com/cse-491/blog/2020/04/20/ranges-transpose.html
#include <numeric>
auto print2D=[](auto rng){for(auto r:rng)print(r);};

TEST_F(RangesTests, ViewsDemo)
{
    namespace rs = std::ranges;
    namespace vs = std::ranges::views;
    std::vector x = {1,2,3,4,5};  // [1,2,3,4,5]
    print(x | vs::drop(2));    // [3,4,5] drop the first two elements of x
    print(x | vs::stride(2));  // [1,3,5] take every other element of x
    print2D(x | vs::chunk(2)); // [[1,2], group elements of x into chunks of length two
                            //  [3,4],
                            //  [5]]
    print(x | vs::chunk(2) | vs::join); // [1,2,3,4,5] join concatenates a range of ranges
    print(x | vs::transform([](auto xi){ return 2*xi; })); // transform maps a lambda over a range [2,4,6,8,10]
    
}

//  a is not a range it is a viewable_range for some reason.  Consequence of chunking?
template <std::ranges::viewable_range V, std::ranges::range R> double my_inner_product(  V&& a,   R&& b, double init)
{
    for (auto&& [ai, bi] : std::views::zip(a, b)) 
        init += ai * bi;   
    return init;
}

TEST_F(RangesTests, Wx)
{
    namespace rs = std::ranges;
    namespace vs = std::ranges::views;
    auto x= vs::iota(1,3+1);
    auto W= vs::iota(1,6+1) | vs::chunk(3); // [1,2,3],[4,5,6]
    auto Wx=W | vs::transform([&](auto row){ return my_inner_product(row, x, 0); });
    print(Wx);
    EXPECT_EQ(Wx[0], 14); // 1*1 + 2*2 + 3*3
    EXPECT_EQ(Wx[1], 32); // 4*1 + 5
   
}

TEST_F(RangesTests, Transpose)
{
    namespace vs = std::ranges::views;
    auto W = vs::iota(1,3*2+1) | vs::chunk(2); // [1,2]
                                            // [3,4]
                                            // [5,6]
    print2D(vs::iota(0,2) | vs::transform([&](int i) { // for each column
    return W | vs::join // concatenate all the rows into a single range
            | vs::drop(i) // remove everything before the 1st element of the ith column
            | vs::stride(2); // take every Nth item to provide a view of the ith column
    })); // [1,3,5]
        // [2,4,6]
}

auto transpose = [](auto rng) {
    namespace rs = std::ranges;
    namespace vs = std::ranges::views;
    auto flat = rng | vs::join;
    int height = rs::distance(rng);
    int width = rs::distance(flat) / height;
    auto inner = [=](int i) mutable {
    return flat | vs::drop(i) | vs::stride(width);
    };
    return vs::iota(0,width) | vs::transform(inner);
};

TEST_F(RangesTests, MatrixMultiplication)
{
    namespace rs = std::ranges;
    namespace vs = std::ranges::views;
    auto X = vs::iota(1,2*3+1) | vs::chunk(3); // [1,2,3]
                                           // [4,5,6]
    auto W = vs::iota(1,3*2+1) | vs::chunk(2); // [1,2]
                                            // [3,4]
                                            // [5,6]
    print2D(X | vs::transform([&](auto xrow) {
    return transpose(W) | vs::transform([=](auto wcol) {
        return my_inner_product(xrow, wcol, 0);
    }); // [22,28]
    }));  // [49,64]
}

