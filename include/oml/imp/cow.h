// File: cow.h  Copy on write, reference counted, array.
#ifndef _cow_h_
#define _cow_h_

// Copyright (1994-2003), Jan N. Reimers

#include "oml/imp/index_t.h"
#include <cassert>



#ifndef OML_USE_STDVEC
template <class T> class cow_array
{
public:
    cow_array(index_t theSize    );
    cow_array(const cow_array& ca);
    ~cow_array(                   );
    cow_array& operator=(const cow_array& ca);

    const T* begin   () const {                         return itsData;}
          T* begin   ()       {if (*itsOwners>1) COW(); return itsData;}
    const T& operator[](index_t i) const {return itsData[i];}
          T& operator[](index_t i)       {if (*itsOwners>1) COW(); return itsData[i];}

    index_t  size() const {return itsSize;}
    int      GetNumOwners() const {return *itsOwners;}

private:
    void Release();
    void COW    ();

          index_t itsSize;
          T*      itsData;
    mutable int*    itsOwners;
};


#ifdef DEBUG
  #define CHECK \
  assert(*itsOwners >0); \
  assert(itsData);       \
  assert(itsSize>=0);
#else
  #define CHECK
#endif
#ifdef WARN_DEEP_COPY
  #include <iostream>
#endif

template <class T> inline cow_array<T>::cow_array(index_t theSize)
  : itsSize  (theSize       )
  , itsData  (new T[itsSize])
  , itsOwners(new int(1)    )
  {
#ifdef WARN_DEEP_COPY
  std::cerr << "*** Copy-on-write array allocating size=" << itsSize << " ***" << std::endl;
#endif
    CHECK
  }

template <class T> inline cow_array<T>::cow_array(const cow_array& ca)
  : itsSize  (ca.itsSize  )
  , itsData  (ca.itsData  )
  , itsOwners(ca.itsOwners)
  {
    assert(itsOwners);
    (*itsOwners)++;
    CHECK
  }

template <class T> inline cow_array<T>::~cow_array()
{
  CHECK
  Release();
}

template <class T> cow_array<T>& cow_array<T>::operator=(const cow_array<T>& ca)
{
  CHECK
  if (this!=&ca)
  {
    Release();
    itsSize  =ca.itsSize;
    itsData  =ca.itsData;
    itsOwners=ca.itsOwners;
    (*itsOwners)++;
  }
  CHECK
  return *this;
}

//
//  Release link to current data.
//
template <class T> inline void cow_array<T>::Release()
{
  CHECK
  if (--(*itsOwners)==0)
  {
    delete [] itsData;
    delete    itsOwners;
#ifdef WARN_DEEP_COPY
  std::cerr << "*** Copy-on-write array de-allocating size=" << itsSize << " ***" << std::endl;
#endif
  }
}

//
//  Deep copy for COW (Copy On Write) symantics.
//
template <class T> void cow_array<T>::COW()
{
  CHECK
#ifdef WARN_DEEP_COPY
  std::cerr << "*** Copy-on-write array doing deep copy, size=" << itsSize << " ***" << std::endl;
#endif

  T* newData=new T[itsSize];
  T* source =itsData;
  T* dest   =newData;
  for (;dest<newData+itsSize; dest++,source++) *dest=*source;
  Release();
  itsData=newData;
  itsOwners=new int(1);
  CHECK
}

#undef CHECK

#endif //USE_STD_VEC

#endif  //_cow_h_
