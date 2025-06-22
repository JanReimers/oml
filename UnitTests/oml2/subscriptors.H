// File: Subscriptors.H  Define dense packing for full,tri,tridiagona,diagonal matrices.
#ifndef _Subscriptors_H_
#define _Subscriptors_H_

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
    typedef std::ranges::iota_view<size_t,size_t> iota_view;
    iota_view nonzero_row_indexes(size_t col) const
    {
        return std::views::iota(size_t(0),nr);
    }
    iota_view nonzero_col_indexes(size_t row) const
    {
        return std::views::iota(size_t(0),nc);
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
    typedef std::ranges::iota_view<size_t,size_t> iota_view;
    iota_view nonzero_row_indexes(size_t col) const
    {
        return std::views::iota(size_t(0),col+1);
    }
    iota_view nonzero_col_indexes(size_t row) const
    {
        return std::views::iota(row,nc);
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
    typedef std::ranges::iota_view<size_t,size_t> iota_view;
    iota_view nonzero_row_indexes(size_t col) const
    {
        return std::views::iota(col,col+1);
    }
    iota_view nonzero_col_indexes(size_t row) const
    {
        return std::views::iota(row,row+1);
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
        assert(nr==nc);
        
        return nr>0 ? 3*nr-2 : 0; // Total number of elements
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
    
    typedef std::ranges::iota_view<size_t,size_t> iota_view;
    iota_view nonzero_row_indexes(size_t col) const
    {
        size_t c0=std::max(size_t(1),col)-1, c1=std::min(nr,col+2);
        return std::views::iota(c0,c1);
    }
    iota_view nonzero_col_indexes(size_t row) const
    {
        size_t r0=std::max(size_t(1),row)-1,r1=std::min(nc,row+2);
        return std::views::iota(r0,r1);
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



#endif //_Subscriptors_H_