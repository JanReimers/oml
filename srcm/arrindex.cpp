module;
#include <cstddef>
#include <iostream>
#include <iomanip>
#include <cmath>
#include <complex>
#include <cassert>
export module oml.ArrIndex;
import oml.VecLimits;
import oml.Shape;
import oml.StreamableObject;
import oml.unop;

export{
//-----------------------------------------------------------------------------
//
//  These macros invoke array index bounds checking if DEBUG is on.
//
#if DEBUG
  #define CHECK(i)\
  assert(i>=0);\
  assert(static_cast<size_t>(i)<size());
#else
  #define CHECK(i)
#endif

template <class T, class Derived,Store M, Shape S> class ArrayIndexable
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

template <class T, class A, Shape S> inline T Sum(const ArrayIndexable<T,A,Full,S>& a)
{
  T ret(0);
  for (const T& ai:a) ret+=ai;
  return ret;
}

template <class T,class A, Store M, Shape S> inline void Fill(ArrayIndexable<T,A,M,S>& arr,T value)
{
  for (T& a:arr) a = value; //Sill fails for complex
}

template <class T,class A, Store M, Shape S> inline void FillLinear(ArrayIndexable<T,A,M,S>& arr,T start, T stop)
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
template <class T,class A, Store M, Shape S> inline A Integrate(const ArrayIndexable<T,A,M,S>& arr,T y0=0)
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

template <class T,class A, Store M, Shape S> inline A Differentiate(const ArrayIndexable<T,A,M,S>& arr)
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
template <class T, class A, Store M, Shape S> inline std::ostream& Write(std::ostream& os,const ArrayIndexable<T,A,M,S>& arr)
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

template <class T, class A, Store M, Shape S> inline std::istream& Read(std::istream& is,ArrayIndexable<T,A,M,S>& arr)
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
template <class T, class Derived,Store M,Shape S> inline \
Derived& Array##NAME (ArrayIndexable<T,Derived,M,S>& a,const T& scalar)\
{\
  for (T& ai:a) ai OP##=scalar;\
  return static_cast<Derived&>(a);\
}\
template <class T,class Derived,Store M,Shape S,class B> inline \
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
template <class T, class A,class B,Store M, Shape S> \
inline bool func(const ArrayIndexable<T,A,M,S>& a, const ArrayIndexable<T,A,M,S>& b) \
{return LogicalII(static_cast<const A&>(a),static_cast<const B&>(b),[](const T& xa,const T&xb){return op;});}\
template <class T, class A,Store M, Shape S> inline bool func(const ArrayIndexable<T,A,M,S>& a, const T& b)\
{return Logical(a,b,[](const T& xa,const T&xb){return op;});}\
template <class T, class A,Store M, Shape S> inline bool func(const T& a, const ArrayIndexable<T,A,M,S>& b)\
{return Logical(a,b,[](const T& xa,const T&xb){return op;});}\

ObScBool(operator==,xa==xb)
ObScBool(operator!=,xa!=xb)

#undef ObScBool
//
//  isnan isinf logical functions
//
inline constexpr bool isnan(const std::complex<double>& c)
{
    return std::isnan(c.real()) || std::isnan(c.imag());
}
inline bool isinf(const std::complex<double>& c)
{
    return std::isinf(c.real()) || std::isinf(c.imag());
}

template <class T, class D, Store M, Shape S> inline bool isnan(const ArrayIndexable<T,D,M,S>& arr)
{
    bool ret=false;
    for (const T& a : arr) if (isnan(a)) {ret=true;break;}
    return ret;
}

template <class T, class D, Store M, Shape S> inline bool isinf(const ArrayIndexable<T,D,M,S>& arr)
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

template <class T, class A, Store M, Shape S> inline T Min(const ArrayIndexable<T,A,M,S>& a)
{
	return MinMax<T,A,OpLT<T>,M,Real,S>::apply(a);
}

template <class T, class A, Store M, Shape S> inline T Max(const ArrayIndexable<T,A,M,S>& a)
{
	return MinMax<T,A,OpGT<T>,M,Real,S>::apply(a);
}

}