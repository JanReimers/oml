#ifndef INDEXABLE_BASE_H_INCLUDED
#define INDEXABLE_BASE_H_INCLUDED

#include "oml/imp/veclimit.h"

template <class Derived, Shape S> class IndexableBase {};

template <class Derived> class IndexableBase<Derived,MatrixShape>
{
    public:
    //
//  Support range based iteration for rows and columns so client code and do
//     for (index_t i:A.rows())
//        for (index_t j:A.cols())
//          {do something with A(i,j)
//
//  For something like a symmtric matrix do
//     for (index_t i:A.rows())
//        for (index_t j:A.cols(i)) //start at i
//          {do something with A(i,j)
//
//
  class index_iterator
  {
    public:
        index_iterator(index_t i) : current{i} {};
        index_iterator operator++(){current++;return (*this);}
        const index_t operator*() const {return current;}
        index_t operator*() {return current;}
        friend bool operator!=(const index_iterator& a, const index_iterator& b) {return a.current!=b.current;}
    private:
        index_t current;
  };

    class iterator_proxy
    {
    public:
        iterator_proxy(const VecLimits& lim) : low(lim.Low), high(lim.High) {};
        iterator_proxy(index_t l, index_t h) : low(l), high(h) {};
        index_iterator begin() const {return low;}
        index_iterator end  () const {return high+1;}
    private:
        index_t low;
        index_t high;
    };

  iterator_proxy rows(         ) const {return iterator_proxy(static_cast<const Derived*>(this)->GetLimits().Row);}
  iterator_proxy cols(         ) const {return iterator_proxy(static_cast<const Derived*>(this)->GetLimits().Col);}
  iterator_proxy rows(index_t j) const {return iterator_proxy(j,static_cast<const Derived*>(this)->GetLimits().Row.High);}
  iterator_proxy cols(index_t i) const {return iterator_proxy(i,static_cast<const Derived*>(this)->GetLimits().Col.High);}
  iterator_proxy array_indices (         ) const {return iterator_proxy(0,static_cast<const Derived*>(this)->size()-1);}
};

template <class Derived> class IndexableBase<Derived,VectorShape>
{
    public:
    //
//  Support range based iteration for rows and columns so client code and do
//     for (index_t i:V)
//          {do something with V(i)
//     for (index_t i:V.all())
//          {do something with V[i]
//
  class index_iterator
  {
    public:
        index_iterator(index_t i) : current{i} {};
        index_iterator operator++(){current++;return (*this);}
        const index_t operator*() const {return current;}
        index_t operator*() {return current;}
        friend bool operator!=(const index_iterator& a, const index_iterator& b) {return a.current!=b.current;}
    private:
        index_t current;
  };

    class iterator_proxy
    {
    public:
        iterator_proxy(const VecLimits& lim) : low(lim.Low), high(lim.High) {};
        iterator_proxy(index_t l, index_t h) : low(l), high(h) {};
        index_iterator begin() const {return low;}
        index_iterator end  () const {return high+1;}
    private:
        index_t low;
        index_t high;
    };

  iterator_proxy indices() const {return iterator_proxy(static_cast<const Derived*>(this)->GetLimits());}
  iterator_proxy array_indices () const {return iterator_proxy(0,static_cast<const Derived*>(this)->size()-1);}
};


#endif // INDEXABLE_BASE_H_INCLUDED
