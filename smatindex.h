// File: smatindex.h  Base class for Glommable Expression Templates of symmetric matricies.
#ifndef _smatindex_h_
#define _smatindex_h_

// Copyright (1994-2003), Jan N. Reimers

#include "oml/indexable.h"

//-------------------------------------------------
//
//  template specialization for symmetruc matricies's.
//
template <class T, class Derived, Data D> class Indexable<T,Derived,Symmetric,D,MatrixShape>
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

  template <class B> Derived& operator+=(const Indexable<T,B,Symmetric,Real,MatrixShape>& b) {return ArrayAdd(*this,b);}
  template <class B> Derived& operator-=(const Indexable<T,B,Symmetric,Real,MatrixShape>& b) {return ArraySub(*this,b);}
  template <class B> Derived& operator+=(const Indexable<T,B,Symmetric,Abstract,MatrixShape>& b) {return MatrixAdd(*this,b);}
  template <class B> Derived& operator-=(const Indexable<T,B,Symmetric,Abstract,MatrixShape>& b) {return MatrixSub(*this,b);}

	class Subscriptor : public Derived::ArraySubscriptor
	{
	 public:
		Subscriptor(Indexable& v) : Derived::ArraySubscriptor(v) {};
	};
	
 protected:
  template <class B> void AssignFrom(const Indexable<T,B,Symmetric,Real    ,MatrixShape>& b) { ArrayAssign(*this,b);}
  template <class B> void AssignFrom(const Indexable<T,B,Symmetric,Abstract,MatrixShape>& b) {MatrixAssign(*this,b);}

 private:
  Indexable& operator=(const Indexable&);
  Indexable(const Indexable&);
};

template <class A> inline bool True(const Indexable<bool,A,Symmetric,Abstract,MatrixShape>& a)                                        
{      
	bool ret(true);
  subsc_t rh=a.GetLimits().Row.High,ch=a.GetLimits().Col.High;
  for (int i=a.GetLimits().Row.Low;i<=rh;i++) 
    for (int j=i;j<=ch;j++)
      ret=ret&&a(i,j);
	return ret;
}                                                                                      

template <class T, class A, Data D> inline T Sum(const Indexable<T,A,Symmetric,D,MatrixShape>& a)                                        
{      
	T ret(0);
  subsc_t rh=a.GetLimits().Row.High,ch=a.GetLimits().Col.High;
  for (int i=a.GetLimits().Row.Low;i<=rh;i++) 
	{
		ret+=a(i,i);
    for (int j=i+1;j<=ch;j++)
      ret+=T(2)*a(i,j);
	}
	return ret;
}                                                                                      

template <class T, class A, Data D> inline 
	std::complex<T> Sum(const Indexable<std::complex<T>,A,Symmetric,D,MatrixShape>& a)                                        
{      
	std::complex<T> ret(0);
  subsc_t rh=a.GetLimits().Row.High,ch=a.GetLimits().Col.High;
  for (int i=a.GetLimits().Row.Low;i<=rh;i++) 
	{
		ret+=a(i,i);
    for (int j=i+1;j<=ch;j++)
      ret+=T(2)*real(a(i,j));
	}
	return ret;
}                                                                                      

template <class T, class A, class Op> class MinMax<T,A,Op,Symmetric,Abstract,MatrixShape>
{
 public:
  static T apply(const Indexable<T,A,Symmetric,Abstract,MatrixShape>& a)
	{
		subsc_t rl=a.GetLimits().Row.Low ,cl=a.GetLimits().Col.Low ;
		subsc_t rh=a.GetLimits().Row.High,ch=a.GetLimits().Col.High;
		T ret=a(rl,cl);
		for (int i=rl;i<=rh;i++) 
			for (int j=i;j<=ch;j++)
		  {
				T ai=a(i,j);
				if (Op::apply(ai,ret)) ret=ai;
			}
		return ret;
	}
};






#endif // _smatindex_h_
