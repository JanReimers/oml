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

//template <class T> inline bool isnan(const Matrix<T>& m)
//{
//    bool ret=false;
//    MatLimits lim=m.GetLimits();
//    for (int i=lim.Row.Low; i<=lim.Row.High; i++)
//        for (int j=lim.Col.Low; j<=lim.Col.High; j++)
//        {
//            if (std::isnan(m(i,j)) || std::isinf(m(i,j)))
//            {
//                ret=true;
//                break;
//            }
//        }
//    return ret;
//}
//
//template <class T> inline bool isnan(const SMatrix<T>& m)
//{
//    bool ret=false;
//    MatLimits lim=m.GetLimits();
//    for (int i=lim.Row.Low; i<=lim.Row.High; i++)
//        for (int j=i; j<=lim.Col.High; j++)
//        {
//            if (std::isnan(m(i,j)) || std::isinf(m(i,j)))
//            {
//                ret=true;
//                break;
//            }
//        }
//    return ret;
//}

//template <class T> inline bool isnan(const Vector<T>& v)
//{
//    bool ret=false;
//    VecLimits lim=v.GetLimits();
//    for (int i=lim.Low; i<=lim.High; i++)
//        if (std::isnan(v(i)) || std::isinf(v(i)))
//        {
//            ret=true;
//            break;
//        }
//    return ret;
//}
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
