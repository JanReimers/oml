// File: smatindex.h  Base class for Glommable Expression Templates of symmetric matricies.
#ifndef _smatindex_h_
#define _smatindex_h_

// Copyright (1994-2003), Jan N. Reimers

#include "oml/imp/matindex.h" //We need to specialization of IndexableBase<D,MatrixShape>

//-------------------------------------------------
//
//  template specialization for symmetruc matricies's.
//
template <class T, class Derived, Data D> class Indexable<T,Derived,Symmetric,D,MatrixShape>
 : public IndexableBase<Derived,MatrixShape>
{
 public:
  explicit Indexable() {};
  ~Indexable() {};

  T  operator()(index_t i,index_t j) const {return static_cast<const Derived*>(this)->operator()(i,j);}

  index_t   size  () const {return static_cast<const Derived*>(this)->size();}
  MatLimits GetLimits() const {return static_cast<const Derived*>(this)->GetLimits();}

 private:
  Indexable& operator=(const Indexable&);
  Indexable(const Indexable&);
};

template <class T, class A, Data D> inline T Sum(const Indexable<T,A,Symmetric,D,MatrixShape>& a)
{
	T ret(0);
  for (index_t i:a.rows())
	{
		ret+=a(i,i);
    for (index_t j:a.cols(i+1))
      ret+=T(2)*a(i,j);
	}
	return ret;
}

template <class T, class A, class Op> class MinMax<T,A,Op,Symmetric,Abstract,MatrixShape>
{
 public:
  static T apply(const Indexable<T,A,Symmetric,Abstract,MatrixShape>& a)
	{
		index_t rl=a.GetLimits().Row.Low ,cl=a.GetLimits().Col.Low ;
		index_t rh=a.GetLimits().Row.High,ch=a.GetLimits().Col.High;
		T ret=a(rl,cl);
		for (index_t i:a.rows())
			for (index_t j:a.cols(i))
		  {
				T ai=a(i,j);
				if (Op::apply(ai,ret)) ret=ai;
			}
		return ret;
	}
};

#endif // _smatindex_h_
