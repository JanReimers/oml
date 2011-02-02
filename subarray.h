// File: subarray.h  Proxy class emulating part of an Array
#ifndef _subarray_h_
#define _subarray_h_

// Copyright (1994-2003), Jan N. Reimers

#include "oml/array.h"

template <class T> class SubArray  
  : public Indexable<T,SubArray<T>,Full,Real,ArrayShape>
  , public Iterable<T,SubArray<T> >
{
 public:
  SubArray(Array<T>& v, index_t size) //assume start=0
    : itsRep(v), itsOffset(0), itsSize(size)
    {
      assert(&itsRep);
      assert(itsSize>=0);
    }

  SubArray(Array<T>& v, subsc_t start, subsc_t stop)
    : itsRep(v), itsOffset(start), itsSize(stop-start+1)
    {
      assert(&itsRep);
      assert(itsOffset>=0);
      assert(itsSize>=0);
    }
  
  SubArray(const SubArray& sv)
    : itsRep(sv.itsRep), itsOffset(sv.itsOffset), itsSize(sv.itsSize) {};

  template <class B> SubArray& operator=(const Indexable<T,B,Full,Real,ArrayShape>& a)
  {
    assert(itsSize==a.size());
    for (int i=0;i<itsSize;i++)
      (*this)[i]=a[i];
    return *this;
  }

  SubArray& operator=(const SubArray& a)
  {
    assert(itsSize==a.size());
    for (int i=0;i<itsSize;i++)
      (*this)[i]=a[i];    
    return *this;
  }

  T  operator[](subsc_t i) const
    {
      assert(i>=0);
      assert(i<itsSize);
      return const_cast<const Array<T>&>(itsRep)[itsOffset+i];
    }
  T& operator[](subsc_t i)
    {
      assert(i>=0);
      assert(i<itsSize);
      return itsRep[itsOffset+i];
    }
  
  index_t size() const {return itsSize;}

  void SetLimits(subsc_t start, subsc_t stop)
  {
    itsOffset=start;
    itsSize=stop-start+1;
    assert(itsOffset>=0);
    assert(itsSize>=0);
  }
#if DEBUG
  #define CHECK(i) assert(i>=0&&i<itsSize)
#else
  #define CHECK(i)
#endif
  class ArraySubscriptor 
  {
   public:
    ArraySubscriptor(Indexable<T,SubArray,Full,Real,ArrayShape>& a) 
      : itsPtr(static_cast<SubArray&>(a).Get()), itsSize(a.size()) {assert(itsPtr);}
    T& operator[](index_t i) {CHECK(i);return itsPtr[i];}
   private:
    T*      itsPtr;
    index_t itsSize;
  };
#undef CHECK


 private:
  friend class Iterable <T,SubArray>;
  friend class ArraySubscriptor;

  const T* Get() const {return &itsRep[0]+itsOffset;} //Required by iterable.
        T* Get()       {return &itsRep[0]+itsOffset;} //Required by iterable.

  Array<T>& itsRep;
  subsc_t   itsOffset;
  subsc_t   itsSize;
};

#endif //_subarray_h_
