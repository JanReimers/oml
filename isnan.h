#ifndef ISNAN_INCLUDED
#define ISNAN_INCLUDED

#include <oml/matrix.h>
#include <oml/smatrix.h>
#include <cmath>
#include <complex>

namespace std
{
template <class T> bool isnan(std::complex<T> c)
{
    return std::isnan(c.real()) || std::isnan(c.imag());
}
template <class T> bool isinf(std::complex<T> c)
{
    return std::isinf(c.real()) || std::isinf(c.imag());
}
}

template <class T, class D> inline bool isnan(const Iterable<T,D>& arr)
{
    bool ret=false;
    typename D::const_iterator i=arr.begin();
    for (;i!=arr.end();i++)
            if (std::isnan(*i) || std::isinf(*i))
        {
            ret=true;
            break;
        }
    return ret;
}
#endif // ISNAN_INCLUDED
