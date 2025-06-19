// File: Ranges.C  Experiment with C++-20 ranges.
//  Following https://gieseanw.wordpress.com/2019/10/20/we-dont-need-no-stinking-expression-templates/

#include "gtest/gtest.h"
#include <ranges>
#include <iostream>
#include <vector>
#include <algorithm>
#include <generator>
#include <valarray>

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
template <std::ranges::range Range> auto inner_product(const Range& a, const Range& b, double init)
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
                            return inner_product(row, x, 0.0);
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
                            return inner_product(row, x, 0.0);
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
            return inner_product(row, x, 0.0);
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
// // x*A = B
// auto operator*(const Vector& x,const Matrix& A)
// {
//     return A | std::views::transform
//     ( [&x](const auto& row)
//         {
//             return inner_product(row, x, 0.0);
//         }
//     );
// }

// auto operator*(const Matrix& A, const Matrix& B)
// {
//     auto product = A | std::views::transform([&B](const auto& row){
//                           return row* B;
//                        });
//     return product;
// }

// TEST_F(RangesTests, MM_Mutiply1)
// {
//     Matrix A = {{1, 2, 3}, {4, 5, 6}};
//     Matrix B = {{7, 8}, {9, 10}, {11, 12}};
//     auto AB = A * B;
//     print(AB[0]);
//     // EXPECT_EQ(AB[0][0],58);
//     // EXPECT_EQ(AB[0][1],64);
//     // EXPECT_EQ(AB[1][0],139);
//     // EXPECT_EQ(AB[1][1],154);
// }