// File: matindex.h  Base class for Glommable Expression Templates.
#ifndef _matindex_h_
#define _matindex_h_

// Copyright (1994-2003), Jan N. Reimers

#include "oml/imp/indexable.h"
#include "oml/imp/memops.h"


//-------------------------------------------------
//
//  template specialization for Matricies's.
//
template <class T, class Derived, Data D> class Indexable<T,Derived,Full,D,MatrixShape>
 : public IndexableBase<Derived,MatrixShape>
{
 public:
  explicit Indexable() {};
  ~Indexable() {};

  T  operator[](index_t n          ) const {return static_cast<const Derived*>(this)->operator[](n);}
  T& operator[](index_t n          )       {return static_cast<      Derived*>(this)->operator[](n);}
  T  operator()(index_t i,index_t j) const {return static_cast<const Derived*>(this)->operator()(i,j);}

  index_t   size  () const {return static_cast<const Derived*>(this)->size();}
  MatLimits GetLimits() const {return static_cast<const Derived*>(this)->GetLimits();}

  Derived& operator+=(T scalar) {return ArrayAdd(*this,scalar);}
  Derived& operator-=(T scalar) {return ArraySub(*this,scalar);}
  Derived& operator*=(T scalar) {return ArrayMul(*this,scalar);}
  Derived& operator/=(T scalar) {return ArrayDiv(*this,scalar);}

  template <class B> Derived& operator+=(const Indexable<T,B,Full,Real,MatrixShape>& b) {return ArrayAdd(*this,b);}
  template <class B> Derived& operator-=(const Indexable<T,B,Full,Real,MatrixShape>& b) {return ArraySub(*this,b);}
  template <class B,Store MB,Data DB> Derived& operator+=(const Indexable<T,B,MB,DB,MatrixShape>& b) {return MatrixAdd(*this,b);}
  template <class B,Store MB,Data DB> Derived& operator-=(const Indexable<T,B,MB,DB,MatrixShape>& b) {return MatrixSub(*this,b);}

	class Subscriptor : public Derived::ArraySubscriptor
	{
	 public:
		Subscriptor(Indexable& m) : Derived::ArraySubscriptor(m) {};
	};

 protected:
  template <class B> void AssignFrom(const Indexable<T,B,Full,Real,MatrixShape>& b) { ArrayAssign(*this,b);}
  template <class B, Store MB, Data DB> void AssignFrom(const Indexable<T,B,MB,DB,MatrixShape>& b) {MatrixAssign(*this,b);}

 private:
  Indexable& operator=(const Indexable&);
  Indexable(const Indexable&);
};
//-------------------------------------------------
//
//  template specialization for abstract Matricies's.
//
template <class T, class Derived> class Indexable<T,Derived,Full,Abstract,MatrixShape>
 : public IndexableBase<Derived,MatrixShape>
{
 public:
  explicit Indexable() {};
  ~Indexable() {};

//  T operator[] Not supported for abstract matricies.
  T operator()(index_t i,index_t j) const {return static_cast<const Derived*>(this)->operator()(i,j);}

  index_t   size  () const {return static_cast<const Derived*>(this)->size();}
  MatLimits GetLimits() const {return static_cast<const Derived*>(this)->GetLimits();}

  Derived& operator+=(T scalar) {return MatrixAdd(*this,scalar);}
  Derived& operator-=(T scalar) {return MatrixSub(*this,scalar);}
  Derived& operator*=(T scalar) {return MatrixMul(*this,scalar);}
  Derived& operator/=(T scalar) {return MatrixDiv(*this,scalar);}

  template <class B,Store MB,Data DB> Derived& operator+=(const Indexable<T,B,MB,DB,MatrixShape>& b) {return MatrixAdd(*this,b);}
  template <class B,Store MB,Data DB> Derived& operator-=(const Indexable<T,B,MB,DB,MatrixShape>& b) {return MatrixSub(*this,b);}

// ArraySubscriptor Not supported for abstract matricies.

 protected:
  template <class B, Store MB, Data DB> void AssignFrom(const Indexable<T,B,MB,DB,MatrixShape>& b) {MatrixAssign(*this,b);}

 private:
  Indexable& operator=(const Indexable&);
  Indexable(const Indexable&);
};

//-------------------------------------------------
//
//  template specialization for diagonal Matricies's.
//
template <class T, class Derived> class Indexable<T,Derived,Diagonal,Abstract,MatrixShape>
 : public IndexableBase<Derived,MatrixShape>
{
 public:
  explicit Indexable() {};
  ~Indexable() {};

//  T operator[] Not supported for abstract matricies.
  T operator()(index_t i,index_t j) const {return static_cast<const Derived*>(this)->operator()(i,j);}

  index_t   size  () const {return static_cast<const Derived*>(this)->size();}
  MatLimits GetLimits() const {return static_cast<const Derived*>(this)->GetLimits();}

  Derived& operator+=(T scalar) {return MatrixAdd(*this,scalar);}
  Derived& operator-=(T scalar) {return MatrixSub(*this,scalar);}
  Derived& operator*=(T scalar) {return MatrixMul(*this,scalar);}
  Derived& operator/=(T scalar) {return MatrixDiv(*this,scalar);}

  template <class B,Store MB,Data DB> Derived& operator+=(const Indexable<T,B,MB,DB,MatrixShape>& b) {return MatrixAdd(*this,b);}
  template <class B,Store MB,Data DB> Derived& operator-=(const Indexable<T,B,MB,DB,MatrixShape>& b) {return MatrixSub(*this,b);}

// ArraySubscriptor Not supported for abstract matricies.

 protected:
  template <class B, Store MB, Data DB> void AssignFrom(const Indexable<T,B,MB,DB,MatrixShape>& b) {MatrixAssign(*this,b);}

 private:
  Indexable& operator=(const Indexable&);
  Indexable(const Indexable&);
};

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

template <class T, class A, class Op, Store M, Data D> class MinMax<T,A,Op,M,D,MatrixShape>
{
public:
    static T apply(const Indexable<T,A,M,D,MatrixShape>& a)
    {
        T ret=a.size()>0 ? a(1,1) : T(0); // Don't try and read a[0] if there is no data in a!
        for (index_t i:a.rows())
        for (index_t j:a.cols())
        {
            T ai=a(i,j);
            if (Op::apply(ai,ret)) ret=ai;
        }
        return ret;
    }
};





#endif // _matindex_h_
