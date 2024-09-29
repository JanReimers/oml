#ifndef TESTS_H_INCLUDED
#define TESTS_H_INCLUDED

#include "gtest/gtest.h"
#include "oml/imp/stream.h"
#include <iostream>
#include <string>

using std::cout;
using std::endl;


template <class Ob> std::string ToString(const Ob& result)
{
    std::stringstream res_stream;
    res_stream << result;
    return res_stream.str();
}

#endif // TESTS_H_INCLUDED
