// File: matindex.h  Base class for Glommable Expression Templates.
#ifndef _matindex_h_
#define _matindex_h_

// Copyright (1994-2003), Jan N. Reimers

#include "oml/imp/indexable.h"

//
// Template specialization provides index iterators for vector shape
//
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
        index_iterator operator++()
        {
            index_iterator ret(*this);
            current++;
            return ret;
        }
        const index_t operator*() const {return current;}
              index_t operator*()       {return current;}
        friend bool operator!=(const index_iterator& a, const index_iterator& b)
        {
            return a.current!=b.current;
        }
    private:
        index_t current;
  };

    class iterator_proxy
    {
    public:
        iterator_proxy(const VecLimits& lim) : low(lim.Low) , high(lim.High)
        {
            if (high<low) high=low-1; //Set to terminate with no iterations
        };
        iterator_proxy(index_t l, index_t h) : low(l), high(h)
        {
            if (high<low) high=low-1; //Set to terminate with no iterations
        };
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

template <class T, class A, class B> class MatrixMMOp;
//-------------------------------------------------
//
//  template specialization for Matrices's.
//
template <class T, class Derived, Data D> class Indexable<T,Derived,Full,D,MatrixShape>
 : public IndexableBase<Derived,MatrixShape>
{
 public:
  explicit Indexable() {};
  ~Indexable() {};

  T  operator()(index_t i,index_t j) const {return static_cast<const Derived*>(this)->operator()(i,j);}

  size_t    size  () const {return static_cast<const Derived*>(this)->size();}
  MatLimits GetLimits() const {return static_cast<const Derived*>(this)->GetLimits();}

 private:
  Indexable& operator=(const Indexable&);
  Indexable(const Indexable&);
};
//-------------------------------------------------
//
//  template specialization for abstract Matrices's.
//
template <class T, class Derived> class Indexable<T,Derived,Full,Abstract,MatrixShape>
 : public IndexableBase<Derived,MatrixShape>
{
 public:
  explicit Indexable() {};
  ~Indexable() {};

  T operator()(index_t i,index_t j) const {return static_cast<const Derived*>(this)->operator()(i,j);}

  size_t    size  () const {return static_cast<const Derived*>(this)->size();}
  MatLimits GetLimits() const {return static_cast<const Derived*>(this)->GetLimits();}

 private:
  Indexable& operator=(const Indexable&);
  Indexable(const Indexable&);
};

//-------------------------------------------------
//
//  template specialization for diagonal Matrices's.
//
template <class T, class Derived> class Indexable<T,Derived,Diagonal,Real,MatrixShape>
 : public IndexableBase<Derived,MatrixShape>
{
 public:
  explicit Indexable() {};
  ~Indexable() {};

  T  operator()(index_t i,index_t j) const {return static_cast<const Derived*>(this)->operator()(i,j);}
  T& operator()(index_t i          )       {return static_cast<      Derived*>(this)->operator()(i);}

  size_t    size  () const {return static_cast<const Derived*>(this)->size();}
  MatLimits GetLimits() const {return static_cast<const Derived*>(this)->GetLimits();}

 private:
  Indexable& operator=(const Indexable&);
  Indexable(const Indexable&);
};

//-------------------------------------------------
//
//  template specialization for abstract Diagonal Matrices's.
//
template <class T, class Derived> class Indexable<T,Derived,Diagonal,Abstract,MatrixShape>
 : public IndexableBase<Derived,MatrixShape>
{
 public:
  explicit Indexable() {};
  ~Indexable() {};

  T operator()(index_t i,index_t j) const {return static_cast<const Derived*>(this)->operator()(i,j);}

  size_t    size  () const {return static_cast<const Derived*>(this)->size();}
  MatLimits GetLimits() const {return static_cast<const Derived*>(this)->GetLimits();}

 private:
  Indexable& operator=(const Indexable&);
  Indexable(const Indexable&);
};

//
//  Create assign functions
//
template <class T, class Derived, Data D, class B, Store MB, Data DB> inline
void MatrixAssign(Indexable<T,Derived,Full,D,MatrixShape>& a,const Indexable<T,B,MB,DB,MatrixShape>& b)
{
#ifdef WARN_DEEP_COPY
  std::cerr << "Doing abstract MatrixAssign n=" << a.size() << std::endl;
#endif
  assert(a.GetLimits()==b.GetLimits());
  typename Derived::Subscriptor s(a);
  #pragma omp parallel for collapse(2)
  for (index_t i=a.GetLimits().Row.Low;i<=a.GetLimits().Row.High;i++)
    for (index_t j=a.GetLimits().Col.Low;j<=a.GetLimits().Col.High;j++)
      s(i,j)=b(i,j);
}

template <class T, class Derived, Data D, class B, Data DB> inline
void MatrixAssign(Indexable<T,Derived,Symmetric,D,MatrixShape>& a,const Indexable<T,B,Symmetric,DB,MatrixShape>& b)
{
#ifdef WARN_DEEP_COPY
  std::cerr << "Doing abstract MatrixAssign n=" << a.size() << std::endl;
#endif
  assert(a.GetLimits()==b.GetLimits());
  typename Derived::Subscriptor s(a);
  #pragma omp parallel for collapse(1)
  for (index_t i=a.GetLimits().Row.Low;i<=a.GetLimits().Row.High;i++)
    for (index_t j=i;j<=a.GetLimits().Col.High;j++)
      s(i,j)=b(i,j);
}

#define OP(NAME,OP) \
template <class T, class Derived,Store M,Data D> inline \
Derived& Matrix##NAME (Indexable<T,Derived,M,D,MatrixShape>& a,const T& scalar)\
{\
  typename Derived::Subscriptor s(a); \
  for (index_t i:a.rows())\
    for (index_t j:a.cols())\
      s(i,j) OP##=scalar;\
  return static_cast<Derived&>(a);\
}\
template <class T, class Derived,Data D> inline \
Derived& Matrix##NAME (Indexable<T,Derived,Symmetric,D,MatrixShape>& a,const T& scalar)\
{\
 typename Derived::Subscriptor s(a); \
 for (index_t i:a.rows())\
    for (index_t j:a.cols())\
      s(i,j) OP##=scalar;\
  return static_cast<Derived&>(a);\
}\
template <class T,class Derived,Store M,Data D,class B,Data DB> inline \
Derived& Matrix##NAME (Indexable<T,Derived,M,D,MatrixShape>& a,\
                  const Indexable<T,B,M,DB,MatrixShape>& b)\
{\
    typename Derived::Subscriptor s(a); \
    for (index_t i:a.rows())\
        for (index_t j:a.cols())\
            s(i,j) OP##=b(i,j);\
  return static_cast<Derived&>(a);\
}\
template <class T,class Derived,Data D,class B,Data DB> inline \
Derived& Matrix##NAME (Indexable<T,Derived,Symmetric,D,MatrixShape>& a,\
                  const Indexable<T,B,Symmetric,DB,MatrixShape>& b)\
{\
  assert(a.GetLimits()==b.GetLimits());\
  typename Derived::Subscriptor s(a); \
  for (index_t i:a.rows())\
    for (index_t j:a.cols(i))\
      s(i,j) OP##=b(i,j);\
  return static_cast<Derived&>(a);\
}

OP(Add,+)
OP(Sub,-)
OP(Mul,*)
OP(Div,/)

#undef OP


template <class T> class Matrix;

template <class T, class A, Store M, Data D> inline
std::ostream& operator<<(std::ostream& os,const Indexable<T,A,M,D,MatrixShape>& a)
{
  return os << Matrix<T>(a);
}

template <class T, class A> inline T Sum(const Indexable<T,A,Full,Abstract,MatrixShape>& a)
{
	T ret(0);
    for (index_t i:a.rows())
        for (index_t j:a.cols())
      ret+=a(i,j);
	return ret;
}

template <class T, class A> inline T Sum(const Indexable<T,A,Diagonal,Abstract,MatrixShape>& a)
{
  T ret(0);
  for (index_t i: a.rows()) ret+=a(i,i);
  return ret;
}


template <class T, class A, class Op, Store M, Data D,Shape S> class MinMax;

template <class T, class A, class Op, Store M> class MinMax<T,A,Op,M,Abstract,MatrixShape>
{
public:
    static T apply(const Indexable<T,A,M,Abstract,MatrixShape>& a)
    {
        int rl=a.GetLimits().Row.Low;
        int cl=a.GetLimits().Col.Low;
        T ret=a.size()>0 ? a(rl,cl) : T(0); // Don't try and read a[0] if there is no data in a!
        for (index_t i:a.rows())
        for (index_t j:a.cols())
        {
            T ai=a(i,j);
            if (Op::apply(ai,ret)) ret=ai;
        }
        return ret;
    }
};

template <class T, class A, Store M> inline T Min(const Indexable<T,A,M,Abstract,MatrixShape>& a)
{
	return MinMax<T,A,OpLT<T>,M,Abstract,MatrixShape>::apply(a);
}

template <class T, class A, Store M> inline T Max(const Indexable<T,A,M,Abstract,MatrixShape>& a)
{
	return MinMax<T,A,OpGT<T>,M,Abstract,MatrixShape>::apply(a);
}

#endif // _matindex_h_
