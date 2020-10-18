// File: matindex.h  Base class for Glommable Expression Templates.
#ifndef _matindex_h_
#define _matindex_h_

// Copyright (1994-2003), Jan N. Reimers

#include "oml/indexable.h"
#include "oml/memops.h"

//-------------------------------------------------
//
//  template specialization for Matricies's.
//
template <class T, class Derived, Data D> class Indexable<T,Derived,Full,D,MatrixShape>
{
 public:
  explicit Indexable() {};
  ~Indexable() {};

  T operator[](index_t n          ) const {return static_cast<const Derived*>(this)->operator[](n);}
  T operator()(subsc_t i,subsc_t j) const {return static_cast<const Derived*>(this)->operator()(i,j);}

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
{
 public:
  explicit Indexable() {};
  ~Indexable() {};

//  T operator[] Not supported for abstract matricies.
  T operator()(subsc_t i,subsc_t j) const {return static_cast<const Derived*>(this)->operator()(i,j);}

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
{
 public:
  explicit Indexable() {};
  ~Indexable() {};

//  T operator[] Not supported for abstract matricies.
  T operator()(subsc_t i,subsc_t j) const {return static_cast<const Derived*>(this)->operator()(i,j);}

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

template <class T> class DMatrix;

template <class T, class A, Store M, Data D> inline
std::ostream& operator<<(std::ostream& os,const Indexable<T,A,M,D,MatrixShape>& a)
{
  return os << DMatrix<T>(a);
}



//----------------------------------------------------------------
//
//  Abstract matrix specializations for some helper functions.
//

template <class A> inline bool True(const Indexable<bool,A,Full,Abstract,MatrixShape>& a)
{
	bool ret(true);
  subsc_t rh=a.GetLimits().Row.High,ch=a.GetLimits().Col.High;
  for (int i=a.GetLimits().Row.Low;i<=rh;i++)
    for (int j=a.GetLimits().Col.Low;j<=ch;j++)
      ret=ret&&a(i,j);
	return ret;
}

template <class T, class A> inline T Sum(const Indexable<T,A,Full,Abstract,MatrixShape>& a)
{
	T ret(0);
  subsc_t rh=a.GetLimits().Row.High,ch=a.GetLimits().Col.High;
  for (int i=a.GetLimits().Row.Low;i<=rh;i++)
    for (int j=a.GetLimits().Col.Low;j<=ch;j++)
      ret+=a(i,j);
	return ret;
}

template <class T, class A, class Op> class MinMax<T,A,Op,Full,Abstract,MatrixShape>
{
 public:
  static T apply(const Indexable<T,A,Full,Abstract,MatrixShape>& a)
	{
		subsc_t rl=a.GetLimits().Row.Low ,cl=a.GetLimits().Col.Low ;
		subsc_t rh=a.GetLimits().Row.High,ch=a.GetLimits().Col.High;
		T ret=a(rl,cl);
		for (int i=rl;i<=rh;i++)
			for (int j=cl;j<=ch;j++)
		  {
				T ai=a(i,j);
				if (Op::apply(ai,ret)) ret=ai;
			}
		return ret;
	}
};






#endif // _matindex_h_
