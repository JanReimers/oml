// File: Vector.H
#ifndef _Vector_H_
#define _Vector_H_

#include <valarray>
#include <ranges>
#include <cassert>

template <std::ranges::viewable_range R> class VectorView
{
public:
    VectorView(R&& r, size_t _n)
        : range(std::forward<R>(r)), n(_n)
    {}

    auto begin()       { return range.begin(); }
    auto end  ()       { return range.end  (); }
    auto begin() const { return range.begin(); }
    auto end  () const { return range.end  (); }
    // many ranges don't support random access with op[].
    size_t size() const { return  n; }

    bool operator==(const std::initializer_list<double>& b) const
    {
        if (size() != b.size()) return false;
        auto it = b.begin();
        for (const auto& v : range) {
            if (v != *it++) return false;
        }
        return true;
    }

    R range; // The underlying range
    size_t n;
};


class Vector
{
public:
    Vector() : data(0) {};
   
    Vector(size_t n)
        : data(n)
    {}

    template <std::ranges::view V> Vector(VectorView<V>& view) //const won't work for views.
        : data(view.size())
    {
        assign_from(view);
    }
    template <std::ranges::view V> Vector(VectorView<V>&& view) //const won't work for views.
        : data(view.size())
    {
        assign_from(view);
    }
    template <std::ranges::view V> Vector& operator=(VectorView<V>&& view) //const won't work for views.
    {   
        if (data.size() == 0)
            data.resize(view.size());
        else
            assert(size()==view.size());

        assign_from(view);
        return *this;
    }

    double operator()(size_t i) const
    {
        assert(i < size());
        return data[i];
    }
    double& operator()(size_t i)
    {
        assert(i < size());
        return data[i];
    }
    size_t size() const { return data.size(); }
    auto  begin()       { return std::begin(data); }
    auto  end  ()       { return std::end  (data); }
    auto  begin() const { return std::begin(data); }
    auto  end  () const { return std::end  (data); }
private:
    template <std::ranges::view V> void assign_from(VectorView<V>& view)
    {
        size_t i=0;
        for (auto& r:view.range) data[i++] = r; //This should be where the lazy evaluation of all the chained views happens.
    }

    std::valarray<double> data;
};



#endif //_Vector_H_