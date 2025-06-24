// File: Matrix.H  Define a matrix class with subscriptor for dense packing.
#ifndef _Matrix_H_
#define _Matrix_H_

#include <oml2/Vector.H>
#include <valarray>
#include <cassert>  

// template <std::ranges::viewable_range R, std::ranges::viewable_range C> class MatrixView
// {
// public:
//     typedef std::ranges::iota_view<size_t,size_t> iota_view;
//     MatrixView(const R& _rows, const C& _cols)
//     : rows(_rows), cols(_cols)
//     {
     
//     }
  
//     size_t size() const { return  nr()*nc(); }
//     size_t nr  () const { return std::ranges::size(rows); }
//     size_t nc  () const { return std::ranges::size(cols); }
// // private:
//     //Maintin two set of ranges one for rows and one for columns.
//     R rows;
//     C cols; 
// };


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
    
    auto row(size_t i) const
    {
        assert(i < subscriptor.nr);
        auto indices=subscriptor.nonzero_col_indexes(i);
        auto v=  indices | std::views::transform([i,this](size_t j){return operator()(i,j);});
        return VectorView<decltype(v)>(std::move(v),indices);
    }
    auto col(size_t j) const
    {
        assert(j < subscriptor.nc);
        auto indices=subscriptor.nonzero_row_indexes(j);
        auto v= indices | std::views::transform([j,this](size_t i){return operator()(i,j);});
        return VectorView<decltype(v)>(std::move(v),indices);
    }
    auto rows() const //Assumes all rows are non-zero.
    {
        return std::views::iota(size_t(0), subscriptor.nr) | std::views::transform([this](size_t i){return row(i);});
    }
    auto cols() const //Assumes all cols are non-zero.
    {
        return std::views::iota(size_t(0), subscriptor.nc) | std::views::transform([this](size_t j){return col(j);});
    }
    auto view () 
    {
        return MatrixView(rows(),cols());
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
// private:
    S subscriptor;
    std::valarray<T> data;
};

#include "oml2/subscriptors.H"

typedef Matrix<double,            FullSubsciptor>            FullMatrix;
typedef Matrix<double, UpperTriangularSubsciptor> UpperTriangularMatrix;
typedef Matrix<double,        DiagonalSubsciptor>        DiagonalMatrix;
typedef Matrix<double,     TriDiagonalSubsciptor>     TriDiagonalMatrix;
// typedef Matrix<double, BandedSubsciptor> BandedMatrix; // 

#endif //_Matrix_H_