// File: Ranges.C  Experiment with C++-20 ranges.
//  Following https://gieseanw.wordpress.com/2019/10/20/we-dont-need-no-stinking-expression-templates/

#include "gtest/gtest.h"
#include <ranges>
#include <iostream>
#include <vector>
#include <algorithm>

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
  
    // cout << p << endl;
    
    // int left = 1;
    // double right = 40.0;
    // auto expr = zip_with(std::plus<>{}, std::ranges::single_view{left}, std::ranges::single_view{right});
    // // print [41]
    // std::cout << expr << std::endl;
        
    // short next_value = 1;
    // auto expr2 = std::ranges::zip_view(std::plus<>{}, expr, std::ranges::single_view{next_value});
        
    // // print [42]
    // std::cout << expr2 << std::endl;
}

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