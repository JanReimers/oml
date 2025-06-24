// File: Matrix.H  Define a matrix class with subscriptor for dense packing.
#ifndef _Matrix_H_
#define _Matrix_H_

#include <oml2/Vector.H>
#include <valarray>
#include <cassert>  

template <typename T, class S> class Matrix
{
    public:
    Matrix(size_t nr, size_t nc) : itsSubscriptor(nc,nr), data(itsSubscriptor.size()) {};
    Matrix(size_t nr, size_t nc, size_t k) : itsSubscriptor(nc,nr,k), data(itsSubscriptor.size()) {};
    Matrix(std::initializer_list<std::initializer_list<T>> init)
        : itsSubscriptor(init.size(), init.begin()->size()), data(itsSubscriptor.size())
    {
        load(init);
    }
    Matrix(std::initializer_list<std::initializer_list<T>> init, size_t k)
    : itsSubscriptor(init.size(), init.begin()->size(),k), data(itsSubscriptor.size())
    {
        load(init);
    }
    T  operator()(size_t i, size_t j) const
    {
        assert(itsSubscriptor.is_stored(i,j));
        return data[itsSubscriptor.offset(i,j)];
    }
    T& operator()(size_t i, size_t j)
    {
        assert(itsSubscriptor.is_stored(i,j));
        return data[itsSubscriptor.offset(i,j)];
    }
    size_t size() const
    {
        return itsSubscriptor.size();
    }
    size_t nr() const { return itsSubscriptor.nr; }
    size_t nc() const { return itsSubscriptor.nc; }
    // auto begin()       { return std::begin(data); }
    // auto end  ()       { return std::end  (data); }
    // auto begin() const { return std::begin(data); }
    // auto end  () const { return std::end  (data); }

    auto row(size_t i) const
    {
        assert(i < itsSubscriptor.nr);
        auto indices=itsSubscriptor.nonzero_col_indexes(i);
        auto v=  indices | std::views::transform([i,this](size_t j){return operator()(i,j);});
        return VectorView<decltype(v)>(std::move(v),indices);
    }
    auto col(size_t j) const
    {
        assert(j < itsSubscriptor.nc);
        auto indices=itsSubscriptor.nonzero_row_indexes(j);
        auto v= indices | std::views::transform([j,this](size_t i){return operator()(i,j);});
        return VectorView<decltype(v)>(std::move(v),indices);
    }
    auto rows() const //Assumes all rows are non-zero.
    {
        return std::views::iota(size_t(0), itsSubscriptor.nr) | std::views::transform([this](size_t i){return row(i);});
    }
    auto cols() const //Assumes all cols are non-zero.
    {
        return std::views::iota(size_t(0), itsSubscriptor.nc) | std::views::transform([this](size_t j){return col(j);});
    }
    S subscriptor() const
    {
        return itsSubscriptor;
    }
    void print() const
    {
        for (size_t i = 0; i < itsSubscriptor.nr; ++i)
        {
            for (size_t j = 0; j < itsSubscriptor.nc; ++j)
            {
                if (itsSubscriptor.is_stored(i, j))
                    std::cout << (*this)(i, j) << " ";
                else
                    std::cout << "0 "; // Print 0 for non-stored elements
            }
            std::cout << std::endl;
        }   
    }
private:
    void load(std::initializer_list<std::initializer_list<T>> init)
    {
        assert(init.size() == itsSubscriptor.nr && "Initializer list size does not match subscriptor row count");
        assert(init.begin()->size() == itsSubscriptor.nc && "Initializer list row size does not match subscriptor column count");
        size_t i = 0;
        for (const auto& row : init)
        {
            size_t j = 0;
            for (const auto& val : row)
            {
                if (itsSubscriptor.is_stored(i, j))
                    data[itsSubscriptor.offset(i, j)] = val;
                else
                    assert(val==0.0);
                ++j;
            }
            ++i;
        }
    }

    S itsSubscriptor;
    std::valarray<T> data;
};

#include "oml2/subscriptors.H"

typedef Matrix<double,            FullSubsciptor>            FullMatrix;
typedef Matrix<double, UpperTriangularSubsciptor> UpperTriangularMatrix;
typedef Matrix<double,        DiagonalSubsciptor>        DiagonalMatrix;
typedef Matrix<double,     TriDiagonalSubsciptor>     TriDiagonalMatrix;
// typedef Matrix<double, BandedSubsciptor> BandedMatrix; // 

#endif //_Matrix_H_