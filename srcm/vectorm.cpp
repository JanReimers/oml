// File: vectorm.cpp  Make a module for Vector<T>
module;
#include <iostream>
#include <iomanip>
#include <cassert>
#include <cstdint>
#include <complex>
#include <vector>

export module oml.vector;
export import oml.StreamableObject;
export import oml.VecLimits;
import oml.Shape;
import oml.unop;
import oml.MixTypes;
import oml.Xpr;
export import oml.Indexable;









// arrindex.h
//-----------------------------------------------------------------------------
//
//  These macros invoke array index bounds checking if DEBUG is on.
//
#if DEBUG
  #include <cassert>
  #define CHECK(i)\
  assert(i>=0);\
  assert(static_cast<size_t>(i)<size());
#else
  #define CHECK(i)
#endif

export template <class T, class Derived,Store M, Shape S> class ArrayIndexable
{
protected: //Can only copy and construct the derived class.
    ArrayIndexable() {};
    ~ArrayIndexable() {};
    ArrayIndexable(const ArrayIndexable&) {};
    ArrayIndexable& operator=(const ArrayIndexable&) {return *this;}

public:
    typedef       T*       iterator;
    typedef const T* const_iterator;

    T operator[](index_t i) const {CHECK(i);return begin()[i];}

    const_iterator begin() const {return static_cast<const Derived*>(this)->priv_begin();}
          iterator begin()       {return static_cast<      Derived*>(this)->priv_begin();}
    const_iterator end  () const {return static_cast<const Derived*>(this)->priv_begin()+size();}
          iterator end  ()       {return static_cast<      Derived*>(this)->priv_begin()+size();}

    size_t  size() const {return static_cast<const Derived*>(this)->size();}

    //  Some overloaded operators
    Derived& operator+=(const T& scalar) {return ArrayAdd(*this,scalar);}
    Derived& operator-=(const T& scalar) {return ArraySub(*this,scalar);}
    Derived& operator*=(const T& scalar) {return ArrayMul(*this,scalar);}
    Derived& operator/=(const T& scalar) {return ArrayDiv(*this,scalar);}

    Derived& operator+=(const ArrayIndexable<T,Derived,M,S>& b) {return ArrayAdd(*this,b);}
    Derived& operator-=(const ArrayIndexable<T,Derived,M,S>& b) {return ArraySub(*this,b);}

    //
    //  Support (index_t i:arr) range iterators over indices
    //
    class index_iterator
    {
    public:
        index_iterator(index_t i) : current{i} {};
        index_iterator operator++(){current++;return (*this);}
        const index_t operator*() const {return current;}
              index_t operator*() {return current;}
        bool operator!=(const index_iterator& b) {return current!=b.current;}
    private:
        index_t current;
    };

    class iterator_proxy
    {
    public:
        iterator_proxy(index_t l, index_t h) : low(l), high(h) {};
        index_iterator begin() const {return low;}
        index_iterator end  () const {return high+1;}
    private:
        index_t low;
        index_t high;
    };

    iterator_proxy arr_indices() const {return iterator_proxy(0,size()-1);}


};

#undef CHECK

export template <class T, class A, Shape S> inline T Sum(const ArrayIndexable<T,A,Full,S>& a)
{
  T ret(0);
  for (const T& ai:a) ret+=ai;
  return ret;
}

export template <class T,class A, Store M, Shape S> inline void Fill(ArrayIndexable<T,A,M,S>& arr,T value)
{
  for (T& a:arr) a = value; //Sill fails for complex
}

export template <class T,class A, Store M, Shape S> inline void FillLinear(ArrayIndexable<T,A,M,S>& arr,T start, T stop)
{
  T del = (stop-start)/(double)(static_cast<A&>(arr).size()-1);
  T val=start;
  for (T& a:arr)
  {
      a = val;
      val+=del;
  }
}

//
//  FE calculus
//
export template <class T,class A, Store M, Shape S> inline A Integrate(const ArrayIndexable<T,A,M,S>& arr,T y0=0)
{
  size_t  n=arr.size();
  A ret(n);
  for (index_t i:arr.indices())
  {
    y0+=arr[i];
    ret[i]=y0;
  }
  return ret;
}

export template <class T,class A, Store M, Shape S> inline A Differentiate(const ArrayIndexable<T,A,M,S>& arr)
{
  index_t n=arr.size();
  A ret(n);
  auto ri=ret.begin();
  *ri=arr[0]; //Save integration constant in case caller needs it.
  for (index_t i=1;i<n;i++,ri++) *ri=arr[i]-arr[i-1];
  return ret;
}

//------------------------------------------------------------------------------
//
//  IO stuff.  Just dumps the data to the stream, ... thats it.
//
export template <class T, class A, Store M, Shape S> inline std::ostream& Write(std::ostream& os,const ArrayIndexable<T,A,M,S>& arr)
{
  assert(os);
  if (StreamableObject::Binary()) for(T b:arr) BinaryWrite(b,os);
  if (StreamableObject::Ascii ()) for(T b:arr) os << b << " ";
  if (StreamableObject::Pretty())
  {
    std::streamsize prec=os.precision();
    std::streamsize wid =os.width();
    os << "{ ";
    for(T b:arr) os << std::setw(wid) << std::setprecision(prec) << b << " ";
    os << "}";
  }
  assert(os);
  return os;
}

export template <class T, class A, Store M, Shape S> inline std::istream& Read(std::istream& is,ArrayIndexable<T,A,M,S>& arr)
{
  assert(is);
  if(StreamableObject::Binary())
    for(T& i:arr) BinaryRead(i,is);
  else
    for(T& i:arr) is >> i;

  assert(is);
  return is;
}
//
//  Create assign functions
//
#define OP(NAME,OP) \
export template <class T, class Derived,Store M,Shape S> inline \
Derived& Array##NAME (ArrayIndexable<T,Derived,M,S>& a,const T& scalar)\
{\
  for (T& ai:a) ai OP##=scalar;\
  return static_cast<Derived&>(a);\
}\
export template <class T,class Derived,Store M,Shape S,class B> inline \
Derived& Array##NAME (ArrayIndexable<T,Derived,M,S>& a, const ArrayIndexable<T,B,M,S>& b)\
{\
    auto bi=b.begin();\
	for (T& ai:a) ai OP##= *bi++;\
	return static_cast<Derived&>(a);\
}\

OP(Add,+)
OP(Sub,-)
OP(Mul,*)
OP(Div,/)

#undef OP


//
//  Logical operators mapped over iterable arrays
//
template <class T, class A, class B, class L, Store M, Shape S>
inline bool LogicalII(const ArrayIndexable<T,A,M,S>& a, const ArrayIndexable<T,B,M,S>& b,const L& lambda)
{
  assert(a.size()==b.size());
  bool ret(true);
  auto bi=b.begin();
  for (const T& ai:a)
  {
    ret = ret && lambda(ai,*bi);
    bi++;
    if (!ret) break;
  }
  return ret;
}

template <class T, class A, class L, Store M, Shape S>
inline bool Logical(const ArrayIndexable<T,A,M,S>& a, const T& b,const L& lambda)
{
  bool ret(true);
  for (const T& i:a) ret = ret && lambda(i,b);
  return ret;
}

template <class T, class A, class L, Store M, Shape S>
inline bool Logical(const T & a, const ArrayIndexable<T,A,M,S>& b,const L& lambda)
{
  bool ret(true);
  for (const T& i:b) ret = ret && lambda(a,i);
  return ret;
}

#define ObScBool(func,op)\
export template <class T, class A,class B,Store M, Shape S> \
inline bool func(const ArrayIndexable<T,A,M,S>& a, const ArrayIndexable<T,A,M,S>& b) \
{return LogicalII(static_cast<const A&>(a),static_cast<const B&>(b),[](const T& xa,const T&xb){return op;});}\
export template <class T, class A,Store M, Shape S> inline bool func(const ArrayIndexable<T,A,M,S>& a, const T& b)\
{return Logical(a,b,[](const T& xa,const T&xb){return op;});}\
export template <class T, class A,Store M, Shape S> inline bool func(const T& a, const ArrayIndexable<T,A,M,S>& b)\
{return Logical(a,b,[](const T& xa,const T&xb){return op;});}\

ObScBool(operator==,xa==xb)
ObScBool(operator!=,xa!=xb)

#undef ObScBool
//
//  isnan isinf logical functions
//
export inline constexpr bool isnan(const std::complex<double>& c)
{
    return std::isnan(c.real()) || std::isnan(c.imag());
}
export inline bool isinf(const std::complex<double>& c)
{
    return std::isinf(c.real()) || std::isinf(c.imag());
}

export template <class T, class D, Store M, Shape S> inline bool isnan(const ArrayIndexable<T,D,M,S>& arr)
{
    bool ret=false;
    for (const T& a : arr) if (isnan(a)) {ret=true;break;}
    return ret;
}

export template <class T, class D, Store M, Shape S> inline bool isinf(const ArrayIndexable<T,D,M,S>& arr)
{
    bool ret=false;
    for (const T& a : arr) if (isinf(a)) {ret=true;break;}
    return ret;
}

//------------------------------------------------------------------
//
//  Max/Min functions.
//  TODO: Try using partial specialization of template functions.
//
template <class T, class A, class Op, Store M, Data D,Shape S> class MinMax;

template <class T, class A, class Op, Store M, Shape S> class MinMax<T,A,Op,M,Real,S>
{
public:
    static T apply(const ArrayIndexable<T,A,M,S>& a)
    {
        T ret=a.size()>0 ? a[0] : T(0); // Don't try and read a[0] if there is no data in a!
        for (index_t i:a.arr_indices())
        {
            T ai=a[i];
            if (Op::apply(ai,ret)) ret=ai;
        }
        return ret;
    }
};

export template <class T, class A, Store M, Shape S> inline T Min(const ArrayIndexable<T,A,M,S>& a)
{
	return MinMax<T,A,OpLT<T>,M,Real,S>::apply(a);
}

export template <class T, class A, Store M, Shape S> inline T Max(const ArrayIndexable<T,A,M,S>& a)
{
	return MinMax<T,A,OpGT<T>,M,Real,S>::apply(a);
}

// vecindex.h
//
// Template specialization provides index iterators for vector shape
//
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
    iterator_proxy indices(index_t i) const 
    {
        return iterator_proxy(i,static_cast<const Derived*>(this)->GetLimits().High);
    }
};

//-------------------------------------------------
//
//  template specialization for Vectors's.
//
template <class T, class Derived, Store M, Data D> class Indexable<T,Derived,M,D,VectorShape>
 : public IndexableBase<Derived,VectorShape>
{
 public:

  T  operator()(index_t n) const {return static_cast<const Derived*>(this)->operator()(n);}
  T& operator()(index_t n)       {return static_cast<      Derived*>(this)->operator()(n);}

  size_t    size  () const {return static_cast<const Derived*>(this)->size();}
  VecLimits GetLimits() const {return static_cast<const Derived*>(this)->GetLimits();}

 protected:
  template <class B> void AssignFrom(const Indexable<T,B,Full,Abstract,VectorShape>& b) {VectorAssign(*this,b);}

  explicit Indexable() {};
  ~Indexable() {};
  Indexable& operator=(const Indexable&) {return *this;}
  Indexable(const Indexable&) {};
};

//--------------------------------------------------------------
//
//  Template specialization for abstract vectors.
//
template <class T, class Derived, Store M> class Indexable<T,Derived,M,Abstract,VectorShape>
 : public IndexableBase<Derived,VectorShape>
{
 public:
  explicit Indexable() {};
  ~Indexable() {};

  T operator()(index_t n) const {return static_cast<const Derived*>(this)->operator()(n);}

  size_t    size  () const {return static_cast<const Derived*>(this)->size();}
  VecLimits GetLimits() const {return static_cast<const Derived*>(this)->GetLimits();}

 private:
  Indexable& operator=(const Indexable&);
  Indexable(const Indexable&);
};

//
//  Create assign functions
//
export template <class T, class Derived,Store M,Data D,class B,Data DB> inline
void VectorAssign(Indexable<T,Derived,M,D,VectorShape>& a,const Indexable<T,B,M,DB,VectorShape>& b)
{
#ifdef WARN_DEEP_COPY
  std::cerr << "Doing abstract VectorAssign n=" << a.size() << std::endl;
#endif
  assert(a.GetLimits()==b.GetLimits());
  typename Derived::Subscriptor s(a);
  for (index_t i:a.indices()) s(i)=b(i);
}

export template <class T, class Derived,Store M,Data D> inline
void VectorAssign(Indexable<T,Derived,M,D,VectorShape>& a,T scalar)
{
#ifdef WARN_DEEP_COPY
  std::cerr << "Doing abstract VectorAssign n=" << a.size() << std::endl;
#endif
  typename Derived::Subscriptor s(a);
  for (index_t i:a.indices()) s(i)=scalar;
}

#define OP(NAME,OP) \
export template <class T, class Derived,Store M,Data D> inline \
Derived& Vector##NAME (Indexable<T,Derived,M,D,VectorShape>& a,const T& scalar)\
{\
  typename Derived::Subscriptor s(a); \
  for (index_t i:a) s(i) OP##=scalar;\
  return static_cast<Derived&>(a);\
}\
export template <class T,class Derived,Store M,Data D,class B,Data DB> inline \
Derived& Vector##NAME (Indexable<T,Derived,M,D,VectorShape>& a,\
                  const Indexable<T,B,M,DB,VectorShape>& b)\
{\
    typename Derived::Subscriptor s(a); \
    for (index_t i:a) s(i) OP##=b(i);\
	return static_cast<Derived&>(a);\
}\

OP(Add,+)
OP(Sub,-)
OP(Mul,*)
OP(Div,/)

#undef OP


//----------------------------------------------------------------
//
//  Abstract vector specializations for some helper functions.
//

export template <class T, class A, Store M> inline
T Sum(const Indexable<T,A,M,Abstract,VectorShape>& a)
{
  T ret(0);
  for (index_t i:a.indices()) ret+=a(i);
  return ret;
}

template <class T, class A, class Op, Store M> class MinMax<T,A,Op,M,Abstract,VectorShape>
{
 public:
  static T apply(const Indexable<T,A,M,Abstract,VectorShape>& a)
  {
    index_t low=a.GetLimits().Low;
    index_t hi =a.GetLimits().High;
    T ret=a(low);
    for (index_t i=low+1;i<=hi;i++)
    {
      T ai=a(i);
      if (Op::apply(ai,ret)) ret=ai;
    }
    return ret;
  }
};

export template <class T, class A, Store M> inline T Min(const Indexable<T,A,M,Abstract,VectorShape>& a)
{
	return MinMax<T,A,OpLT<T>,M,Abstract,VectorShape>::apply(a);
}

export template <class T, class A, Store M> inline T Max(const Indexable<T,A,M,Abstract,VectorShape>& a)
{
	return MinMax<T,A,OpGT<T>,M,Abstract,VectorShape>::apply(a);
}

// minmax.h
export template <class T> inline       T& Max(      T& a,       T& b)  {return a > b ? a : b;}
export template <class T> inline const T& Max(const T& a, const T& b)  {return a > b ? a : b;}
export template <class T> inline       T& Min(      T& a,       T& b)  {return a < b ? a : b;}
export template <class T> inline const T& Min(const T& a, const T& b)  {return a < b ? a : b;}

// tstram.h
export template <class A> class TStreamableObject;

export template <class A> std::ostream& operator<<(std::ostream& os, const TStreamableObject<A>& o);
export template <class A> std::istream& operator>>(std::istream& os,       TStreamableObject<A>& o);

// Allows template based streaming.  IO methods for the parent type A are called directly.
//  ***No virtual dispatch***!  If you need virtual dispatch see the PMStreamableObject class.
export template <class A> class TStreamableObject
: public StreamableObject
{
 public:
  friend std::ostream& operator<< <>(std::ostream& os, const TStreamableObject& o);
  friend std::istream& operator>> <>(std::istream& os,       TStreamableObject& o);
//  friend std::ostream& operator<<(std::ostream& os, const TStreamableObject* o) {return os << *o;}
//  friend std::istream& operator>>(std::istream& is,       TStreamableObject* o) {return is >> *o;}
  static A* Factory(std::istream& is);
};


export template <class A> inline std::ostream& operator<<(std::ostream& os, const TStreamableObject<A>& o)
{
  o.WriteHeader(os,typeid(A).name());
  return static_cast<const A&>(o).Write(os);
}

export template <class A> inline std::istream& operator>>(std::istream& is,TStreamableObject<A>& o)
{
  StreamableObject::Mode current=o.ReadHeader(is,typeid(A).name());
  static_cast<A&>(o).Read (is);
  o.SetOutputMode(current); //Restore to previous state.
  return is;
}

template <class A> inline A* TStreamableObject<A>::Factory(std::istream& is)
{
  //  CheckName(is,typeid(A).name()); Read operation should do this.
  A* ret=new A;
  is >> ret;
  return ret;
}


// cow.h
#ifndef OML_USE_STDVEC
export template <class T> class cow_array
{
public:
    cow_array(size_t theSize     );
    cow_array(const cow_array& ca);
    ~cow_array(                   );
    cow_array& operator=(const cow_array& ca);

    const T* begin   () const {                         return itsData;}
          T* begin   ()       {if (*itsOwners>1) COW(); return itsData;}
    const T& operator[](index_t i) const {return itsData[i];}
          T& operator[](index_t i)       {if (*itsOwners>1) COW(); return itsData[i];}

    size_t   size() const {return itsSize;}
    int      GetNumOwners() const {return *itsOwners;}

private:
    void Release();
    void COW    ();

          size_t  itsSize;
          T*      itsData;
    mutable int*  itsOwners;
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

template <class T> inline cow_array<T>::cow_array(size_t theSize)
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

// ran250.h
#include <complex>
#include <cassert>
#include <climits>
#ifdef HAVE_STDINT_H
  #include <stdint.h>
#else
  #ifdef HAVE_INTTYPES_H
    #include <inttypes.h>
  #else
    #ifdef  HAVE_SYS_TYPES_H
      #include <sys/types.h>
      typedef  __uint32_t uint32_t;
    #endif
  #endif
#endif

template <class T> union Int2RealUnion;

template <> union Int2RealUnion<float>
{		   	// used to access floats as unsigneds
  float rep;
  uint32_t u;
};

template <> union Int2RealUnion<double>
{		   	// used to access doubles as unsigneds
  double rep;
  uint32_t u[2];
};

template <class T> class Int2Real;

template <> class Int2Real<float>
{
 public:
  Int2Real();
  float convert(uint32_t l) const;
 private:
  Int2RealUnion<float> itsMantissa;
};

template <> class Int2Real<double>
{
 public:
  Int2Real();
  double convert(uint32_t l1,uint32_t l2) const;
 private:
  Int2RealUnion<double> itsMantissa;
};


inline float Int2Real<float>::convert(uint32_t l) const
{
  Int2RealUnion<float> result;
  result.rep = 1.0;
  result.u   |= (l & itsMantissa.u);
  result.rep -= 1.0;
  assert( result.rep < 1.0 && result.rep >= 0);
  return result.rep;
}

inline double Int2Real<double>::convert(uint32_t l1,uint32_t l2) const
{
  Int2RealUnion<double> result;
  result.rep = 1.0;
  result.u[0] |= (l1 & itsMantissa.u[0]);
  result.u[1] |= (l1 & itsMantissa.u[1]);
  result.rep -= 1.0;
  assert( result.rep < 1.0 && result.rep >= 0);
  return result.rep;
}

class TwoTap
{
  unsigned int   Next,Tap1,Tap2,Mask;
  long*          itsArray;

  Int2Real<float > floatConverter;
  Int2Real<double> doubleConverter;

public:
  TwoTap(unsigned int tap1,unsigned int tap2);
 ~TwoTap();
  const char* Name();

  long GetNext()
  {
    ++Next;
    return itsArray[Next&Mask]=
//      (itsArray[(Next-Tap1)&Mask]%(SHRT_MAX-2))
//       *
//      (itsArray[(Next-Tap2)&Mask]%(SHRT_MAX-2));
      itsArray[(Next-Tap1)&Mask]^itsArray[(Next-Tap2)&Mask];
  }

  float  GetNextFloat () {return  floatConverter.convert(GetNext());}
  double GetNextDouble() {return doubleConverter.convert(GetNext(),GetNext());}


};


class FourTap
{
  unsigned int   Next,Tap1,Tap2,Tap3,Tap4,Mask;
  long*          itsArray;

  Int2Real<float > floatConverter;
  Int2Real<double> doubleConverter;


public:
  FourTap(unsigned int tap1,unsigned int tap2,unsigned int tap3,unsigned int tap4);
 ~FourTap();
  const char* Name();

  long GetNext()
  {
    ++Next;
    return itsArray[Next&Mask]=
      itsArray[(Next-Tap1)&Mask]^
      itsArray[(Next-Tap2)&Mask]^
      itsArray[(Next-Tap3)&Mask]^
      itsArray[(Next-Tap4)&Mask];
  }
  float  GetNextFloat () {return  floatConverter.convert(GetNext());}
  double GetNextDouble() {return doubleConverter.convert(GetNext(),GetNext());}
};

#ifdef HAVE_STDLIB_H
#include <stdlib.h>
#endif
#include <time.h>
#ifdef HAVE_LIMITS_H
#include <limits.h>
#endif
#include <stdio.h>
#include <cmath>
#include <cassert>
#ifdef HAVE_SYS_TYPES_H
#include <sys/types.h>
#endif

//
//	Additive number generator. This method is presented in Volume II
//	of The Art of Computer Programming by Knuth. I've coded the algorithm
//	and have added the extensions by Andres Nowatzyk of CMU to randomize
//	the result of algorithm M a bit	by using an LCG & a spatial
//	permutation table.
//
//	The version presented uses the same constants for the LCG that Andres
//	uses (chosen by trial & error). The spatial permutation table is
//	the same size (it's based on word size). This is for 32-bit words.
//
//	The ``auxillary table'' used by the LCG table varies in size, and
//	is chosen to be the the smallest power of two which is larger than
//	twice the size of the state table.
//

class ACG
{
  uint32_t initialSeed;	// used to reset generator
  int initialTableEntry;

  uint32_t *state;
  uint32_t *auxState;
  short stateSize;
  short auxSize;
  uint32_t lcgRecurr;
  short j;
  short k;
public:
  ACG(uint32_t seed = 0, int size = 55);
  ~ACG();
  uint32_t asLong();
  void reset();
};

double INorm=1.0/LONG_MAX;
double LNorm=1.0/LONG_MAX;

TwoTap::TwoTap(unsigned int tap1,unsigned int tap2)
  : Next(0)
  , Tap1(tap1)
  , Tap2(tap2)
  , Mask(0)
  , itsArray(0)
{
  unsigned int max=Max(Tap1,Tap2);
  unsigned int n=(unsigned int)pow(2.0,(int)floor(log(max)/log(2.0)+1.0));
  itsArray=new long[n];
  Mask=n-1;
  ACG Rand(time(0)); //Use this guy to boot strap R(250,103).
  for (unsigned int i=0;i<n;i++) itsArray[i]=Rand.asLong();
}

TwoTap::~TwoTap()
{
  delete [] itsArray;
}

const char* TwoTap::Name()
{
    std::ostringstream os;
    os << "R(" << Tap1 << "," << Tap2 << ",*)";
    return os.str().c_str();
}

FourTap::FourTap(unsigned int tap1,unsigned int tap2,unsigned int tap3,unsigned int tap4)
  : Next(0)
  , Tap1(tap1)
  , Tap2(tap2)
  , Tap3(tap3)
  , Tap4(tap4)
  , Mask(0)
  , itsArray(0)
{
  unsigned int max=Max(Max(Tap1,Tap2),Max(Tap3,Tap4));
  unsigned int n=(unsigned int)pow(2.0,(int)floor(log(max)/log(2.0)+1.0));
  itsArray=new long[n];
  Mask=n-1;
  ACG Rand(time(0)); //Use this guy to boot strap R(250,103).
  for (unsigned int i=0;i<n;i++) itsArray[i]=Rand.asLong();
}

FourTap::~FourTap()
{
  delete [] itsArray;
}

const char* FourTap::Name()
{
    std::ostringstream os;
    os << "R(" << Tap1 << "," << Tap2 << "," << Tap3 << "," << Tap4 << ")";
    return os.str().c_str();
}


FourTap GlobalRandomNumberGenerator(471,1586,6988,9869);
//TwoTap GlobalRandomNumberGenerator(1063,1279);
//TwoTap GlobalRandomNumberGenerator(103,250);

//
//	The following is a hack that I attribute to
//	Andres Nowatzyk at CMU. The intent of the loop
//	is to form the smallest number 0 <= x < 1.0,
//	which is then used as a mask for two longwords.
//	this gives us a fast way way to produce double
//	precision numbers from longwords.
//
//	I know that this works for IEEE and VAX floating
//	point representations.
//
//	A further complication is that gnu C will blow
//	the following loop, unless compiled with -ffloat-store,
//	because it uses extended representations for some of
//	of the comparisons. Thus, we have the following hack.
//	If we could specify #pragma optimize, we wouldn't need this.
//
Int2Real<double>::Int2Real()
{
  assert (sizeof(double) == 2 * sizeof(uint32_t));
  Int2RealUnion<double> t;
#if _IEEE == 1

  t.rep = 1.5;
  if ( t.u[1] == 0 )
  {		// sun word order?
    t.u[0] = 0x3fffffff;
    t.u[1] = 0xffffffff;
  }
  else
  {
    t.u[0] = 0xffffffff;	// encore word order?
    t.u[1] = 0x3fffffff;
  }

#else
  volatile double x = 1.0; // volatile needed when fp hardware used,
  // and has greater precision than memory doubles
  double y = 0.5;
  do
  {			    // find largest fp-number < 2.0
    t.rep = x;
    x += y;
    y *= 0.5;
  }
  while (x != t.rep && x < 2.0);

#endif
  // set doubleMantissa to 1 for each doubleMantissa bit
  itsMantissa.rep = 1.0;
  itsMantissa.u[0] ^= t.u[0];
  itsMantissa.u[1] ^= t.u[1];
}

Int2Real<float>::Int2Real()
{
  assert (sizeof(float) == sizeof(uint32_t));
  Int2RealUnion<float> t;
#if _IEEE == 1
  t.u = 0x3fffffff;
#else
  t.rep=0.0; //avoid valgrind warning.
  volatile float x = 1.0; // volatile needed when fp hardware used,
  // and has greater precision than memory floats
  float y = 0.5;
  do
  {			    // find largest fp-number < 2.0
    t.rep = x;
    x += y;
    y *= 0.5;
  }
  while (x != t.rep && x < 2.0);
#endif
  // set singleMantissa to 1 for each singleMantissa bit
  itsMantissa.rep = 1.0;
  itsMantissa.u ^= t.u;

}




extern double  INorm;
extern double  LNorm;

export {
template <class T> T      OMLRand();
template <class T> T      OMLRandPos();
template <class T> double OMLRandScale(T max);

template <> inline int    OMLRand<int>   () {return GlobalRandomNumberGenerator.GetNext();}
template <> inline long   OMLRand<long>  () {return GlobalRandomNumberGenerator.GetNext();}
template <> inline float  OMLRand<float> ()
{
  return GlobalRandomNumberGenerator.GetNextFloat();
}
template <> inline double OMLRand<double>()
{
  return GlobalRandomNumberGenerator.GetNextDouble();
}

template <> inline std::complex<double> OMLRand<std::complex<double> >()
{
  return std::complex<double>(OMLRand<double>(),OMLRand<double>());
}

template <> inline int    OMLRandPos<int>   () {return OMLRand<int> ()&0x7fffffff;}
template <> inline long   OMLRandPos<long>  () {return OMLRand<long>()&0x7fffffff;}
template <> inline float  OMLRandPos<float> () {return OMLRand<float>();}
template <> inline double OMLRandPos<double>() {return OMLRand<double>();}
template <> inline std::complex<double> OMLRandPos<std::complex<double> >()
{
  return std::complex<double>(OMLRandPos<double>(),OMLRandPos<double>());
}


template <> inline double OMLRandScale<int>   (int    max) {return INorm*max;}
template <> inline double OMLRandScale<long>  (long   max) {return INorm*max;}
template <> inline double OMLRandScale<float> (float  max) {return max;}
template <> inline double OMLRandScale<double>(double max) {return max;}
template <> inline double OMLRandScale<std::complex<double> >(std::complex<double> max) {return real(max);}

// random.h
/*! \file random.h
  \brief Routines for filling OML containers with random numbers.

  Random number are supplied by a 2 tap xor random number generator. Very fast!
 */
//! Fill with random numbers, range = (-1,1).
template <class T,class A,Store M, Shape S> void FillRandom(ArrayIndexable<T,A,M,S>& arr)
{
  for (T& i:arr) i=OMLRand<T>();
}

//! Fill with positive random numbers, range = [0,1).
template <class T,class A,Store M, Shape S> void FillRandomPositive(ArrayIndexable<T,A,M,S>& arr)
{
  for (T& i:arr) i=OMLRandPos<T>();
}

//! Fill with positive random numbers, range = [-max,max).
template <class T,class A,Store M, Shape S> void FillRandom(ArrayIndexable<T,A,M,S>& arr,T max)
{
  double scale=OMLRandScale<T>(max);
  for (T& i:arr) i = (T)(OMLRand<T>() * scale);
}

//! Fill with positive random numbers, range = [0,max).
template <class T,class A,Store M, Shape S> void FillRandomPositive(ArrayIndexable<T,A,M,S>& arr,T max)
{
  double scale=OMLRandScale<T>(max);
  for (T& i:arr) i = (T)(OMLRandPos<T>() * scale);
}
} //export block
// vector3d.h
const double pi=acos(-1.0);

//-----------------------------------------------------------------------------
/*! \class Vector3D vector3d.h oml/vector3d.h
  \brief Very light weight 3D %Vector class with lots of overloaded operators.

  Structure for real space three vectors.  Most algeabraic operators have
  been overloaded.  So Vector3Ds should behave like any intrinsic data type.
  Lots of overloaded operators form 3 element vectors.
  - +, -, +=, -=, ==, !=
  - \a a \c % \a b is the vector cross product.
  - \c * is a dot product.
  - /, *=, /= for scalers.
  - <, >, <=, >= all compare magnitudes.
  - Unary + and -
  - \a a| \a b returns the angle in radians between \a a and \a b.
  - \a a|| \a b returns the angle in degrees between \a a and \a b.
  - ! \a a returns the magnitude of \a a.
  - ~ \a a returns a unit vector in the direction of \a a.
  - For complex values \c conj, \c real, \c imag and \c norm are defined.

  This class is much more efficient that \c Vector<double>(3) would be.
  You will need to include \c io3d.h to get \c op<< and \c op>> for IO.
  \nosubgrouping
*/
export template <class T> class Vector3D
{
 public:
  /*! \name Constructors/Assignment*/
  //@{

  //! Default contructor, \c x=y=z=0.
  Vector3D(                 ): x(  0),y(  0),z(  0) {};
  //! Contruct from individual components.
  Vector3D(T _x,T _y,T _z   ): x( _x),y( _y),z( _z) {};
  //! Copy constructor.
  Vector3D(const Vector3D& v): x(v.x),y(v.y),z(v.z) {};
  //! Construct from another data type.
  template <class T1> Vector3D(T1 _x,T1 _y,T1 _z    ) : x( _x),y( _y),z( _z) {}
  //! Construct form another Vector3D type.
  template <class T1> Vector3D(const Vector3D<T1>& v) : x(v.x),y(v.y),z(v.z) {}

  //! Assign
  Vector3D& operator =(const Vector3D& v) {x=v.x;y=v.y;z=v.z;return *this;}
  //! Assign form another Vector3D type.
  template <class T1> Vector3D& operator =(const Vector3D<T1>& v) {x=v.x;y=v.y;z=v.z;return *this;}
  //@}

 ~Vector3D() {};

  //! Element access
  const T& operator()(index_t i) const {return (&x)[i-1];}
  //! Element access
  T& operator()(index_t i)       {return (&x)[i-1];}

 


  /*! \name Coordinates*/
  //@{
  T x; //!< \a x coordinate.
  T y; //!< \a y coordinate.
  T z; //!< \a z coordinate.
  //@}
};

export {
//-----------------------------------------------------------------------------
//
//  Binary algeabra.
//


template <class T1, class T2> inline
auto operator +(const Vector3D<T1>& a,const Vector3D<T2>& b)
{
  typedef typename ReturnType<T1,T2>::RetType TR;
  return Vector3D<TR>(a.x+b.x,a.y+b.y,a.z+b.z);
}

template <class T1, class T2> inline
auto operator -(const Vector3D<T1>& a,const Vector3D<T2>& b)
{
  typedef typename ReturnType<T1,T2>::RetType TR;
  return Vector3D<TR>(a.x-b.x,a.y-b.y,a.z-b.z);
}

template <class T1, class T2> inline
auto operator *(const Vector3D<T1>& a,const Vector3D<T2>& b)
{
  //typedef typename ReturnType<T1,T2>::RetType TR;
  return (a.x*b.x+a.y*b.y+a.z*b.z);
}

template <class T> inline
Vector3D<T> operator *(const Vector3D<T>& a,const T F)
{
  return Vector3D<T>(a.x*F,a.y*F,a.z*F);
}

template <class T> inline
Vector3D<T> operator *(const T F,const Vector3D<T>& a)
{
  return Vector3D<T>(a.x*F,a.y*F,a.z*F);
}

template <class T> inline
Vector3D<T> operator /(const Vector3D<T>& a,const T F)
{
  return Vector3D<T>(a.x/F,a.y/F,a.z/F);
}

template <class T1, class T2> inline
auto operator %(const Vector3D<T1>& a,const Vector3D<T2>& b)
{
  typedef typename ReturnType<T1,T2>::RetType TR;
  return Vector3D<TR>( a.x%b.x, a.y%b.y, a.z%b.z );
}

template <class T1, class T2> inline
auto Cross(const Vector3D<T1>& a,const Vector3D<T2>& b)
{
  typedef typename ReturnType<T1,T2>::RetType TR;
  return Vector3D<TR>( a.y*b.z-a.z*b.y, a.z*b.x-a.x*b.z, a.x*b.y-a.y*b.x );
}

//
//  A operator= B overloads for binary operators
//

template <class T1, class T2> inline
Vector3D<T1>& operator +=(Vector3D<T1>& a,const Vector3D<T2>& b)
{
  a.x+=b.x; a.y+=b.y; a.z+=b.z;
  return a;
}

template <class T1, class T2> inline
Vector3D<T1>& operator -=(Vector3D<T1>& a,const Vector3D<T2>& b)
{
  a.x-=b.x; a.y-=b.y; a.z-=b.z;
  return a;
}

template <class T> inline
Vector3D<T>& operator *=(Vector3D<T>& a,const T F)
{
  a.x*=F;a.y*=F;a.z*=F;
  return a;
}

template <class T> inline
Vector3D<T>& operator /=(Vector3D<T>& a,const T F)
{
  a.x/=F;a.y/=F;a.z/=F;
  return a;
}

//-----------------------------------------------------------------------------
//
//  Relational operators.
//
template <class T> inline
bool operator ==(const Vector3D<T>& a,const Vector3D<T>& b)
{
  return ((a.x==b.x)&&(a.y==b.y)&&(a.z==b.z));
}

template <class T> inline
bool operator !=(const Vector3D<T>& a,const Vector3D<T>& b)
{
  return !(a==b);
}

template <class T> inline
bool operator > (const Vector3D<T>& a,const Vector3D<T>& b)
{
  return (norm(a) > norm(b));
}

template <class T> inline
bool operator < (const Vector3D<T>& a,const Vector3D<T>& b)
{
  return (norm(a) < norm(b));
}

template <class T> inline
bool operator >=(const Vector3D<T>& a,const Vector3D<T>& b)
{
  return (norm(a) >= norm(b));
}

template <class T> inline
bool operator <=(const Vector3D<T>& a,const Vector3D<T>& b)
{
  return (norm(a) <= norm(b));
}

//-----------------------------------------------------------------------------
//
//  Unary operators
//
template <class T> inline
Vector3D<T>  operator -(const Vector3D<T>& a)
{
  return Vector3D<T>(-a.x,-a.y,-a.z);
}

template <class T> inline
Vector3D<T>  operator +(const Vector3D<T>& a)
{
  return a;
}

//-----------------------------------------------------------------------------
//
//  Angle between two vectors (Radians).
//
template <class T> inline
T angle(const Vector3D<T>& a,const Vector3D<T>& b)
{
  return (T)acos( normalize(a) * normalize(b) );
}

//-----------------------------------------------------------------------------
//
//  Angle between two vectors (Degrees).
//
template <class T> inline
T angle_degrees(const Vector3D<T>& a,const Vector3D<T>& b)
{
  return static_cast<T>(acos( normalize(a) * normalize(b) )/pi*180.0);
}

//-----------------------------------------------------------------------------
//
//  Magnitude and normalize.
//
template <class T> inline
T norm(const Vector3D<T>& a)
{
  return (T)sqrt(a*a);
}

template <class T> inline
Vector3D<T>  normalize(const Vector3D<T>& a)  //normalize.
{
  return a/norm(a);
}

//--------------------------------------------------------------------
//
//  Specialized templates for complex data types.
//
template <class T> inline Vector3D<std::complex<T> > conj(const Vector3D<std::complex<T> >& v)
{
  return Vector3D<std::complex<T> >(conj(v.x),conj(v.y),conj(v.z));
}

template <class T> inline Vector3D<T> real(const Vector3D<std::complex<T> >& v)
{
  return Vector3D<T>(real(v.x),real(v.y),real(v.z));
}


template <class T> inline Vector3D<T> imag(const Vector3D<std::complex<T> >& v)
{
  return Vector3D<T>(imag(v.x),imag(v.y),imag(v.z));
}

} //export block



// vector.h
/*! \class Vector vector.h oml/vector.h
  \brief Numerical container with FORTRAN array symatics.

  Vectors have element indexes ranging from 1...n (by default) and are indexed using the
  (i) syntax, just like FORTRAN arrays.  The base index can also be changed just like in FORTRAN.
  Vectors support efficient copy on write (COW) symmantices, just like Arrays.

  \b Math:

  The special operators for Vector are
  - \c operator* does a dot product (not a direct multiply).
  - for complex value \c V*V will do V-dot-conj(V).
  - use \c DirectMultiply(V,V) to get a Vector of element by element products.
  - \c operator/ is is only allowed for scalars.
  - \c Magnitude(V) and \c !V both return the magnitude of V.
  - \c Normalize(V) returns a unit vector with the same direction as V.

  \nosubgrouping
*/

export enum class FillType {Zero,Random,Unit};

export template <class T> class Vector
  : public Indexable<T,Vector<T>,Full,Real,VectorShape>
  , public ArrayIndexable<T,Vector<T>,Full,VectorShape>
  , public TStreamableObject<Vector<T> >
{
 public:
  typedef Indexable<T,Vector<T>,Full,Real,VectorShape> IndexableT;
  typedef ArrayIndexable <T,Vector<T>,Full     ,VectorShape> IterableT;
  typedef Ref<T,IndexableT,VectorShape> RefT;
  typedef Ref<T,IterableT ,VectorShape> RefAT;

  /*! \name Constructors/Assignment
  Copy constructor and op=(Vector) are automaically supplied by the compiler.
  */
  //@{
  //!< Vector with size=0;
           Vector() : Vector<T>(VecLimits(1,0)) {};
  //!<  All elements are un-initiated, low index is 1.
  explicit Vector(size_t  size) : Vector<T>(VecLimits(1,size)) {};
  explicit Vector(size_t  size, FillType ft) : Vector<T>(VecLimits(1,size),ft) {};
  explicit Vector(size_t  size, const T&  fillValue) : Vector<T>(VecLimits(1,size),fillValue) {};
  //!<  Specify lower and upper index.
//  explicit Vector(index_t l,index_t h) : Vector<T>(VecLimits(l,h)) {};
  explicit Vector(index_t l,index_t h, FillType ft) : Vector<T>(VecLimits(l,h),ft) {};
  explicit Vector(index_t l,index_t h, const T&  fillValue) : Vector<T>(VecLimits(l,h),fillValue) {};
  Vector(const VecLimits&); //!<  Specify lower and upper index.
  Vector(const VecLimits&, FillType ft); //!<  Specify lower and upper index.
  Vector(const VecLimits&, const T&  fillValue       ); //!<  Specify lower and upper index.
  //! Allows construction from an expression template.
  template <class B,Data D> Vector(const Indexable<T,B,Full,D,VectorShape>&);
  //! Allows assignment from an expression template.
  template <class B,Data D> Vector& operator=(const Indexable<T,B,Full,D,VectorShape>&);

  Vector(const Vector& m);
  Vector& operator=(const Vector&);

#ifdef OML_MOVE_OPS
  Vector(Vector&& m);
  Vector& operator=(Vector&&);
#endif

  VecLimits ReBase(int low);
  VecLimits ReBase(const VecLimits& lim);
  //@}
  void Fill(FillType);
  void Fill(const T& fillValue);
  void FillRandom();

  std::ostream& Write(std::ostream&) const;
  std::istream& Read (std::istream&)      ;
  // Need to disambiguate from expression version.
  friend std::ostream& operator<<(std::ostream& os,const Vector& a)
  {
    return os << static_cast<const TStreamableObject<Vector<T> >& >(a);
  }


  /*! \name Subscripting operators
    FORTRAN style 1-D array subscripting.
    If DEBUG is defined, every index will be checked that it is in range.
    For fast write access with op() make a subscriptor, then the COW check is only done once
    during construction of the Subscriptor.
  */
  //@{
  //! const element acces operator, fast and \e cannot trigger a COW operation.
  T  operator()(index_t) const;
  //! non-const version can trigger a COW operation, and checks for this with every access.
  T& operator()(index_t)      ;
  //@}

  size_t    size     () const; //!<Returns number elements in the Vector.
  VecLimits GetLimits() const; //!<Returns lower and upper limits in a structure.
  index_t   GetLow   () const; //!<Returns lower index.
  index_t   GetHigh  () const; //!<Returns upper index.

  /*! \name Resize functions
    The data can be optionally preserved.
  */
  //@{
  void SetLimits(const VecLimits&, bool preserve=false); //!<Resize from new limits.
  void SetLimits(size_t          , bool preserve=false); //!<Resize from new size.
  void SetLimits(index_t,index_t , bool preserve=false); //!<Resize from new limits.
  //@}

  //! Does V(i)=V(index[i]) for i=low...high.  Used for sorting.
  void   ReIndex(const std::vector<index_t>&);

  /*! \name Sub-vector functions  */
  //@{
  Vector SubVector(const VecLimits&) const; //!< From limits.
  Vector SubVector(size_t          ) const; //!< First n elements.
  Vector SubVector(index_t,index_t ) const; //!< From limits.
  //@}

  /*! \name Iterators.
    Iterators should be STL compatable.
   */
  //@{
  //! Read only iterator.
  typedef typename IterableT::const_iterator  const_iterator;
  //! Read/write iterator.
  typedef typename IterableT::iterator iterator;
  //@}

#if DEBUG
  #define CHECK(i) assert(itsLimits.CheckIndex(i))
#else
  #define CHECK(i)
#endif
  class Subscriptor
  {
   public:
    Subscriptor(Indexable<T,Vector,Full,Real,VectorShape>& a)
      : itsLimits(a.GetLimits())
      , itsPtr(static_cast<Vector*>(&a)->priv_begin()-itsLimits.Low)
      {assert(itsPtr);}

    T& operator()(index_t i) {CHECK(i);return itsPtr[i];}

   private:
    VecLimits itsLimits;
    T*        itsPtr;
  };
#undef CHECK

 private:
  friend class      Indexable<T,Vector,Full,Real,VectorShape>;
  friend class ArrayIndexable<T,Vector,Full     ,VectorShape>;
  friend class Subscriptor;

  const T* priv_begin() const {return &*itsData.begin();} //Required by iterable.
        T* priv_begin()       {return &*itsData.begin();} //Required by iterable.
  void  Check () const; //Check internal consistency between limits and cow.

  VecLimits    itsLimits; //Manages the upper and lower vector limits.
#ifdef OML_USE_STDVEC
  std::vector<T> itsData;
#else
  cow_array<T> itsData;   //Copy-On-Write array for the data.
#endif
};

export {
//-------------------------------------------------------------------------
//
//  Some special operators only for vectors.
//
template <class T, class A, Store M, Data D> inline
std::ostream& operator<<(std::ostream& os,const Indexable<T,A,M,D,VectorShape>& a)
{
  return os << Vector<T>(a);
}


// The define op* as a vector dot (inner) product.
template <class TA, class TB, class A, class B, Store MA, Store MB, Data DA, Data DB> inline
auto operator*(const Indexable<TA,A,MA,DA,VectorShape>& a, const Indexable<TB,B,MB,DB,VectorShape>& b)
{
    typedef typename ReturnType<TA,TB>::RetType TR;
    TR ret(0);
    for (index_t k:a.indices()) ret+=a(k)*b(k);
    return ret;
}

// Calculate the magnitue of a vector.
template <class T, class A, Store M, Data D> inline
T Magnitude(const Indexable<T,A,M,D,VectorShape>& a)
{
  return sqrt(a*a);
}

// Overload op! for magnitue.
template <class T, class A, Store M, Data D> inline
T operator!(const Indexable<T,A,M,D,VectorShape>& a)
{
  return Magnitude(a);
}

// Rescale so that v*v==1.0 .
template <class T, class A, Store M, Data D> inline
void Normalize(Indexable<T,A,M,D,VectorShape>& v)
{
  v/=!v;
}
} //export block
//-----------------------------------------------------------------------------
//
//  Macro, which expands to an index checking function call,
//  when DEBUG is on
//
#if DEBUG
  #define CHECK(i) itsLimits.CheckIndex(i)
#else
  #define CHECK(i)
#endif

template <class T> inline T Vector<T>::operator()(index_t i) const
{
  CHECK(i);
  return itsData[itsLimits.Offset(i)];
}

template <class T> inline T& Vector<T>::operator()(index_t i)
{
  CHECK(i);
  return itsData[itsLimits.Offset(i)];
}

#undef CHECK

#ifdef DEBUG
#define CHECK Check()
#else
#define CHECK
#endif

template <class T> inline Vector<T>::Vector(const VecLimits& lim,FillType ft)
  : Vector<T>(lim)
  {
    Fill(ft);
  }

template <class T> inline Vector<T>::Vector(const VecLimits& lim,const T& fillValue)
  : Vector<T>(lim)
  {
     Fill(fillValue);
  }

//
//  All other constructors should delegate to this one.
//
template <class T> inline Vector<T>::Vector(const VecLimits& lim)
  : itsLimits(lim          )
  , itsData  (lim.size())
  {
    CHECK;
  }

template <class T> inline Vector<T>::Vector(const Vector<T>& v)
  : itsLimits(v.itsLimits)
  , itsData  (v.itsData)
  {}

template <class T> inline Vector<T>& Vector<T>::operator=(const Vector<T>& v)
{
  itsLimits=v.itsLimits;
  itsData  =v.itsData;
//  std::cout << "Vector<T> move op=" << std::endl;
  return *this;
}
#ifdef OML_MOVE_OPS

template <class T> inline Vector<T>::Vector(Vector<T>&& v)
  : itsLimits(std::move(v.itsLimits))
  , itsData  (std::move(v.itsData))
  {
//    std::cout << "Vector<T> move constructor m.itsData.size()=" << m.itsData.size() << std::endl;
  }

template <class T> inline Vector<T>& Vector<T>::operator=(Vector<T>&& v)
{
  itsLimits=std::move(v.itsLimits);
  itsData  =std::move(v.itsData);
//  std::cout << "Vector<T> move op=" << std::endl;
  return *this;
}
#endif

template <class T> inline VecLimits Vector<T>::ReBase(int low)
{
    VecLimits oldLimits=itsLimits;
    itsLimits.ReBase(low);
    return oldLimits;
}

template <class T> inline VecLimits Vector<T>::ReBase(const VecLimits& lim)
{
    VecLimits oldLimits=itsLimits;
    itsLimits.ReBase(lim.Low);
    return oldLimits;
}

template <class T> inline void Vector<T>::Fill(FillType ft)
{
    switch (ft)
    {
        case (FillType::Zero)   : Fill(T(0));break;
        case (FillType::Random) : ::FillRandom(*this);break;
        case (FillType::Unit)   : Fill(T(1));break;
    }
}

template <> inline Vector3D<double> OMLRand<Vector3D<double> >()
{
  return Vector3D<double>(OMLRand<double>(),OMLRand<double>(),OMLRand<double>());
}

template <>  inline void Vector<Vector3D<double> >::Fill(FillType ft)
{
    switch (ft)
    {
        case (FillType::Zero)   : Fill(Vector3D<double>(double(0),double(0),double(0)));break;
        case (FillType::Random) : ::FillRandom(*this);break;
        case (FillType::Unit)   : Fill(Vector3D<double>(double(1),double(1),double(1)));break;
    }
}

template <class T> inline void Vector<T>::Fill(const T& fillValue)
{
    ::Fill(*this,fillValue);
}

template <class T> inline void Vector<T>::FillRandom()
{
    ::FillRandom(*this);
}


template <class T> inline size_t  Vector<T>::size() const
{
  return GetLimits().size();
}



template <class T> template <class B,Data D> inline
Vector<T>::Vector(const Indexable<T,B,Full,D,VectorShape>& v)
  : itsLimits(v.GetLimits())
  , itsData  (itsLimits.size())
  {
    this->AssignFrom(v); //Choose op[i] or op(i) depending on whether v is abstract.
    CHECK;
  }

#undef CHECK

template <class T> inline  VecLimits Vector<T>::GetLimits() const
{
  return itsLimits;
}

template <class T> inline index_t Vector<T>::GetLow() const
{
  return itsLimits.Low;
}

template <class T> inline index_t Vector<T>::GetHigh() const
{
  return itsLimits.High;
}

template <class T> inline void Vector<T>::SetLimits(size_t size,bool preserve)
{
  SetLimits(VecLimits(size),preserve);
}

template <class T> inline void Vector<T>::SetLimits(index_t l,index_t h,bool preserve)
{
  SetLimits(VecLimits(l,h),preserve);
}

template <class T> inline Vector<T> Vector<T>::SubVector(size_t size) const
{
  return SubVector(VecLimits(size));
}

template <class T> inline Vector<T> Vector<T>::SubVector(index_t l,index_t h) const
{
  return SubVector(VecLimits(l,h));
}

template <class T> template <class B, Data D> inline
Vector<T>& Vector<T>::operator=(const Indexable<T,B,Full,D,VectorShape>& v)
{
  if (size()==0) SetLimits(v.GetLimits());
  this->AssignFrom(v); //Choose op[i] or op(i) depending on whether v is abstract.
  return *this;
}

#ifdef DEBUG
#define CHECK Check()
#else
#define CHECK
#endif



//-----------------------------------------------------------------------------
//
//  Changing size and/or limits.
//
template <class T> void Vector<T>::SetLimits(const VecLimits& theLimits, bool preserve)
{
  theLimits.Check();
  if (itsLimits!=theLimits)
  {
    if (preserve)
    {
      Vector<T> dest(theLimits);

      index_t low =Max(GetLow() ,theLimits.Low );   //Limits of old and new
      index_t high=Min(GetHigh(),theLimits.High);   //data overlap.

      Subscriptor source(*this);         //Destination subscriptor.

      for (index_t i=low;i<=high;i++) dest(i)=source(i); //Transfer any overlaping data.
      *this=dest;
    }
    else
    {
      *this=Vector<T>(theLimits);
    }
  }
  CHECK;
}

//  Concatenate two Vectors into a new Vector.
export template <class T> Vector<T> operator&(const Vector<T>& a, const Vector<T>& b)
{
  typedef typename Vector<T>::const_iterator CI;
  VecLimits newlim(a.GetLimits().Low,a.GetLimits().Low+a.size()+b.size()-1);
  Vector<T> ret(newlim);
  typename Vector<T>::iterator i=ret.begin();
  for (CI ab=a.begin();ab!=a.end();ab++,i++) *i=*ab;
  for (CI bb=b.begin();bb!=a.end();bb++,i++) *i=*bb;
  return ret;
}

export template <class T> void Vector<T>::ReIndex(const std::vector<index_t>& index)
{
  assert(size()==index.size());

  std::vector<index_t>::const_iterator b=index.begin();
  Vector<T>                      dest(GetLimits());
  iterator                       i=dest.begin();
  for (;b!=index.end();b++,i++) *i=(*this)[*b];
  *this=dest;
}

export template <class T> Vector<T> Vector<T>::SubVector(const VecLimits& lim) const
{
  assert(GetLimits().Low  <= lim.Low );
  assert(GetLimits().High >= lim.High);
  Vector<T> ret(lim);
  Subscriptor r(ret);
  for (index_t i=ret.GetLimits().Low;i<=ret.GetLimits().High;i++) r(i)=(*this)(i);
  return ret;
}


#if DEBUG
//-----------------------------------------------------------------------------
//
//  Internal consistency check.
//
template <class T> void Vector<T>::Check() const
{
  assert(itsLimits.Check());
  assert(itsLimits.size()==itsData.size());
}
#endif //DEBUG
//-----------------------------------------------------------------------------
//
//  IO
//
template <class T> std::ostream& Vector<T>::Write(std::ostream& os) const
{
  assert(os);
  if (!this->Pretty())
    os << GetLimits();
   else
  {
    std::streamsize prec=os.precision();
    std::streamsize wid =os.width();
    os << GetLimits();
    os << std::setw(wid) << std::setprecision(prec);
  }
  return ::Write(os,*this);
}

template <class T> std::istream& Vector<T>::Read(std::istream& is)
{
  assert(is);
  VecLimits lim;
  is >> lim;
  if (size()==0)
    SetLimits(lim);
   else
    assert(GetLimits()==lim);

  ::Read(is,*this);
  CHECK;
  assert(is);
  return is;
}



#undef CHECK

// template class Vector<double>;
// template class Vector<std::complex<double>>;



//
//	This is an extension of the older implementation of Algorithm M
//	which I previously supplied. The main difference between this
//	version and the old code are:
//
//		+ Andres searched high & low for good constants for
//		  the LCG.
//
//		+ theres more bit chopping going on.
//
//	The following contains his comments.
//
//	agn@UNH.CS.CMU.EDU sez..
//
//	The generator below is based on 2 well known
//	methods: Linear Congruential (LCGs) and Additive
//	Congruential generators (ACGs).
//
//	The LCG produces the longest possible sequence
//	of 32 bit random numbers, each being unique in
//	that sequence (it has only 32 bits of state).
//	It suffers from 2 problems: a) Independence
//	isnt great, that is the (n+1)th number is
//	somewhat related to the preceding one, unlike
//	flipping a coin where knowing the past outcomes
//	dont help to predict the next result.  b)
//	Taking parts of a LCG generated number can be
//	quite non-random: for example, looking at only
//	the least significant byte gives a permuted
//	8-bit counter (that has a period length of only
//	256).  The advantage of an LCA is that it is
//	perfectly uniform when run for the entire period
//	length (and very uniform for smaller sequences
//	too, if the parameters are chosen carefully).
//
//	ACGs have extremly long period lengths and
//	provide good independence.  Unfortunately,
//	uniformity isnt not too great. Furthermore, I
//	didnt find any theoretically analysis of ACGs
//	that addresses uniformity.
//
//	The RNG given below will return numbers
//	generated by an LCA that are permuted under
//	control of a ACG. 2 permutations take place: the
//	4 bytes of one LCG generated number are
//	subjected to one of 16 permutations selected by
//	4 bits of the ACG. The permutation a such that
//	byte of the result may come from each byte of
//	the LCG number. This effectively destroys the
//	structure within a word. Finally, the sequence
//	of such numbers is permuted within a range of
//	256 numbers. This greatly improves independence.
//
//
//  Algorithm M as describes in Knuths "Art of Computer Programming",
//	Vol 2. 1969
//  is used with a linear congruential generator (to get a good uniform
//  distribution) that is permuted with a Fibonacci additive congruential
//  generator to get good independence.
//
//  Bit, byte, and word distributions were extensively tested and pass
//  Chi-squared test near perfect scores (>7E8 numbers tested, Uniformity
//  assumption holds with probability > 0.999)
//
//  Run-up tests for on 7E8 numbers confirm independence with
//  probability > 0.97.
//
//  Plotting random points in 2d reveals no apparent structure.
//
//  Autocorrelation on sequences of 5E5 numbers (A(i) = SUM X(n)*X(n-i),
//	i=1..512)
//  results in no obvious structure (A(i) ~ const).
//
//  Except for speed and memory requirements, this generator outperforms
//  random() for all tests. (random() scored rather low on uniformity tests,
//  while independence test differences were less dramatic).
//
//  AGN would like to..
//  thanks to M.Mauldin, H.Walker, J.Saxe and M.Molloy for inspiration & help.
//
//  And I would (DGC) would like to thank Donald Kunth for AGN for letting me
//  use his extensions in this implementation.
//

//
//	Part of the table on page 28 of Knuth, vol II. This allows us
//	to adjust the size of the table at the expense of shorter sequences.
//

static int randomStateTable[][3] = {
{3,7,16}, {4,9, 32}, {3,10, 32}, {1,11, 32}, {1,15,64}, {3,17,128},
{7,18,128}, {3,20,128}, {2,21, 128}, {1,22, 128}, {5,23, 128}, {3,25, 128},
{2,29, 128}, {3,31, 128}, {13,33, 256}, {2,35, 256}, {11,36, 256},
{14,39,256}, {3,41,256}, {9,49,256}, {3,52,256}, {24,55,256}, {7,57, 256},
{19,58,256}, {38,89,512}, {17,95,512}, {6,97,512}, {11,98,512}, {-1,-1,-1} };

//
// spatial permutation table
//	RANDOM_PERM_SIZE must be a power of two
//

#define RANDOM_PERM_SIZE 64
uint32_t randomPermutations[RANDOM_PERM_SIZE] = {
0xffffffff, 0x00000000,  0x00000000,  0x00000000,  // 3210
0x0000ffff, 0x00ff0000,  0x00000000,  0xff000000,  // 2310
0xff0000ff, 0x0000ff00,  0x00000000,  0x00ff0000,  // 3120
0x00ff00ff, 0x00000000,  0xff00ff00,  0x00000000,  // 1230

0xffff0000, 0x000000ff,  0x00000000,  0x0000ff00,  // 3201
0x00000000, 0x00ff00ff,  0x00000000,  0xff00ff00,  // 2301
0xff000000, 0x00000000,  0x000000ff,  0x00ffff00,  // 3102
0x00000000, 0x00000000,  0x00000000,  0xffffffff,  // 2103

0xff00ff00, 0x00000000,  0x00ff00ff,  0x00000000,  // 3012
0x0000ff00, 0x00000000,  0x00ff0000,  0xff0000ff,  // 2013
0x00000000, 0x00000000,  0xffffffff,  0x00000000,  // 1032
0x00000000, 0x0000ff00,  0xffff0000,  0x000000ff,  // 1023

0x00000000, 0xffffffff,  0x00000000,  0x00000000,  // 0321
0x00ffff00, 0xff000000,  0x00000000,  0x000000ff,  // 0213
0x00000000, 0xff000000,  0x0000ffff,  0x00ff0000,  // 0132
0x00000000, 0xff00ff00,  0x00000000,  0x00ff00ff   // 0123
};

//
//	SEED_TABLE_SIZE must be a power of 2
//
#define SEED_TABLE_SIZE 32
static uint32_t seedTable[SEED_TABLE_SIZE] = {
0xbdcc47e5, 0x54aea45d, 0xec0df859, 0xda84637b,
0xc8c6cb4f, 0x35574b01, 0x28260b7d, 0x0d07fdbf,
0x9faaeeb0, 0x613dd169, 0x5ce2d818, 0x85b9e706,
0xab2469db, 0xda02b0dc, 0x45c60d6e, 0xffe49d10,
0x7224fea3, 0xf9684fc9, 0xfc7ee074, 0x326ce92a,
0x366d13b5, 0x17aaa731, 0xeb83a675, 0x7781cb32,
0x4ec7c92d, 0x7f187521, 0x2cf346b4, 0xad13310f,
0xb89cff2b, 0x12164de1, 0xa865168d, 0x32b56cdf
};

//
//	The LCG used to scramble the ACG
//
//
// LC-parameter selection follows recommendations in
// "Handbook of Mathematical Functions" by Abramowitz & Stegun 10th, edi.
//
// LC_A = 251^2, ~= sqrt(2^32) = 66049
// LC_C = result of a long trial & error series = 3907864577
//

static const uint32_t LC_A = 66049;
static const uint32_t LC_C = 3907864577U;
static inline uint32_t LCG(uint32_t x)
{
    return( x * LC_A + LC_C );
}


ACG::ACG(uint32_t seed, int size)
{
    int l;
    initialSeed = seed;

    //
    //	Determine the size of the state table
    //

    for (l = 0;
	 randomStateTable[l][0] != -1 && randomStateTable[l][1] < size;
	 l++);

    if (randomStateTable[l][1] == -1) {
	l--;
    }

    initialTableEntry = l;

    stateSize = randomStateTable[ initialTableEntry ][ 1 ];
    auxSize = randomStateTable[ initialTableEntry ][ 2 ];

    //
    //	Allocate the state table & the auxillary table in a single malloc
    //

    state = new uint32_t[stateSize + auxSize];
    auxState = &state[stateSize];

    reset();
}

//
//	Initialize the state
//
void
ACG::reset()
{
    uint32_t u;

    if (initialSeed < SEED_TABLE_SIZE) {
	u = seedTable[ initialSeed ];
    } else {
	u = initialSeed ^ seedTable[ initialSeed & (SEED_TABLE_SIZE-1) ];
    }


    j = randomStateTable[ initialTableEntry ][ 0 ] - 1;
    k = randomStateTable[ initialTableEntry ][ 1 ] - 1;

    int i;
    for(i = 0; i < stateSize; i++) {
	state[i] = u = LCG(u);
    }

    for (i = 0; i < auxSize; i++) {
	auxState[i] = u = LCG(u);
    }

    k = u % stateSize;
    int tailBehind = (stateSize - randomStateTable[ initialTableEntry ][ 0 ]);
    j = k - tailBehind;
    if (j < 0) {
	j += stateSize;
    }

    lcgRecurr = u;

    assert(sizeof(double) == 2 * sizeof(int32_t));
}

ACG::~ACG()
{
    if (state) delete [] state;
    state = 0;
    // don't delete auxState, it's really an alias for state.
}

//
//	Returns 32 bits of random information.
//

uint32_t
ACG::asLong()
{
    uint32_t result = state[k] + state[j];
    state[k] = result;
    j = (j <= 0) ? (stateSize-1) : (j-1);
    k = (k <= 0) ? (stateSize-1) : (k-1);

    short int auxIndex = (result >> 24) & (auxSize - 1);
    uint32_t auxACG = auxState[auxIndex];
    auxState[auxIndex] = lcgRecurr = LCG(lcgRecurr);

    //
    // 3c is a magic number. We are doing four masks here, so we
    // do not want to run off the end of the permutation table.
    // This insures that we have always got four entries left.
    //
    uint32_t *perm = & randomPermutations[result & 0x3c];

    result =  *(perm++) & auxACG;
    result |= *(perm++) & ((auxACG << 24)
			   | ((auxACG >> 8)& 0xffffff));
    result |= *(perm++) & ((auxACG << 16)
			   | ((auxACG >> 16) & 0xffff));
    result |= *(perm++) & ((auxACG <<  8)
			   | ((auxACG >> 24) &   0xff));

    return(result);
}


