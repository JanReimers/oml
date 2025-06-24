// File: Vector.H
#ifndef _Vector_H_
#define _Vector_H_

#include <valarray>
#include <ranges>
#include <cassert>

template <std::ranges::viewable_range R> class VectorView
{
public:
    typedef std::ranges::iota_view<size_t,size_t> iota_view;
   
    VectorView(R&& r, const iota_view& indices)
    : range(std::forward<R>(r)), itsIndices(indices)
    {
        assert(size() == std::ranges::size(range) && "VectorView stop index out of range");
    }
    VectorView(R& r, const iota_view& indices)
    : range(std::forward<R>(r)), itsIndices(indices)
    {
        assert(size() == std::ranges::size(range) && "VectorView stop index out of range");
    }
     VectorView(R& r)
    : range(std::forward<R>(r)), itsIndices(size_t(0),range.size())
    {
        assert(size() == std::ranges::size(range) && "VectorView stop index out of range");
    }
    VectorView(const R& r, const iota_view& indices)
    : range(r), itsIndices(indices)
    {
        assert(size() == std::ranges::size(range) && "VectorView stop index out of range");
    }

    auto begin()       { return std::ranges::begin(range); }
    auto end  ()       { return std::ranges::end  (range); }
    auto begin() const { return std::ranges::begin(range); }
    auto end  () const { return std::ranges::end  (range); }
    // many ranges don't support random access with op[].
   
    size_t size() const { return  itsIndices.size(); }
    iota_view indices() const { return itsIndices; }

    bool operator==(const std::initializer_list<double>& b) const
    {
        if (size() != b.size()) return false;
        auto it = b.begin();
        for (const auto& v : range) {
            if (v != *it++) return false;
        }
        return true;
    }
private:
    R range; // only includes data for the non-zero portion of the vector.
    iota_view itsIndices;
    
};


template <class T> class Vector
{
public:
    Vector() : data(0) {};
   
    Vector(size_t n)
        : data(n)
    {}
    Vector(std::initializer_list<T> init)
        : data(init.size())
    {
        size_t i = 0;
        for (const auto& val : init) {
            data[i++] = val;
        }
    }

    template <std::ranges::view V> Vector(const VectorView<V>& view) //const won't work for views.
        : data(view.size())
    {
        assign_from(view);
    }
    // template <std::ranges::view V> Vector(VectorView<V>&& view) //const won't work for views.
    //     : data(view.size())
    // {
    //     assign_from(view);
    // }

    template <std::ranges::view V> Vector& operator=(const VectorView<V>&& view) //const won't work for views.
    {   
        if (data.size() == 0)
            data.resize(view.size());
        else
            assert(size()==view.size());

        assign_from(view);
        return *this;
    }

    T operator()(size_t i) const
    {
        assert(i < size());
        return data[i];
    }
    T& operator()(size_t i)
    {
        assert(i < size());
        return data[i];
    }
    size_t size() const { return data.size(); }
    auto  begin()       { return std::begin(data); }
    auto  end  ()       { return std::end  (data); }
    auto  begin() const { return std::begin(data); }
    auto  end  () const { return std::end  (data); }
    typedef std::ranges::iota_view<size_t,size_t> iota_view;
    iota_view indices() const { return iota_view(size_t(0), size()); } //Full view of indices

    auto view() const
    {
        return VectorView<decltype(data)>(data,iota_view(size_t(0), size())); //Full view
    }
    auto view(const iota_view& indices) const //possibly a partial view
    {
        auto v =indices | std::views::transform([this](size_t i){ return data[i]; });
        return VectorView<decltype(v)>(std::move(v), indices);
    }
    template <std::ranges::viewable_range R> operator VectorView<R>() const
    {
        return VectorView<R>(data, iota_view(size_t(0), size()));
    }
private:
    template <std::ranges::view V> void assign_from(const VectorView<V>& view)
    {
        size_t i=0;
        for (auto r:view) data[i++] = r; //This should be where the lazy evaluation of all the chained views happens.
    }

    std::valarray<T> data;
};



#endif //_Vector_H_