// File: subvector.h  Proxy class emulating part of a vector

// Copyright (1994-2003), Jan N. Reimers

#include "oml/vector.h"

template <class T> class SubVector
  : public Indexable<T,SubVector<T>,Full,Abstract,VectorShape>
//  , public Iterable<T,SubVector<T> >
//  , public TStreamableObject<SubVector<T> >
{
 public:
  SubVector(Vector<T>& v, index_t start, index_t size)
    : itsRep(v), itsOffset(0), itsSize(size), itsLow(start)
    {
      assert(&itsRep);
      assert(start>=itsRep.GetLow());
      assert(start+itsSize<itsRep.GetHigh());
    }

  SubVector(Vector<T>& v, index_t start, index_t size,index_t low)
    : itsRep(v), itsOffset(start-low), itsSize(size), itsLow(low)
    {
      assert(&itsRep);
      assert(start>=itsRep.GetLow());
      assert(start+itsSize-1<=itsRep.GetHigh());
    }

  SubVector(const SubVector& sv)
    : itsRep(sv.itsRep), itsOffset(sv.itsOffset), itsSize(sv.itsSize), itsLow(sv.itsLow) {};

  template <class B> SubVector& operator=(const Indexable<T,B,Full,Real,VectorShape>& a)
  {
    assert(itsSize==a.size());
    for (index_t i=0;i<itsSize;i++)
      (*this)(i+itsLow)=a(i+itsLow);
    return *this;
  }

  SubVector& operator=(const SubVector& a)
  {
    assert(itsSize==a.size());
    for (index_t i=0;i<itsSize;i++)
      (*this)(i+itsLow)=a(i+itsLow);
    return *this;
  }

  T  operator()(index_t i) const
    {
      assert(i>=itsLow);
      assert(i<itsLow+itsSize);
      return const_cast<const Vector<T>&>(itsRep)(itsOffset+i);
    }
  T& operator()(index_t i)
    {
      assert(i>=itsLow);
      assert(i<itsLow+itsSize);
      return itsRep(itsOffset+i);
    }

  index_t size() const {return itsSize;}

 private:
  const T* Get() const {return itsRep.Get()+itsOffset;} //Required by iterable.
        T* Get()       {return itsRep.Get()+itsOffset;} //Required by iterable.

  Vector<T>& itsRep;
  index_t    itsOffset;
  index_t    itsSize;
  index_t    itsLow;  //Lower index
};
