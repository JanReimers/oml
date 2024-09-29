// File: Vector.h  Interface base class for std::vector objects.
#ifndef _Vector_h_
#define _Vector_h_

// Copyright (1994-2005), Jan N. Reimers

#include "oml/imp/veclimit.h"
#include "oml/imp/arrindex.h"
#include "oml/imp/vecindex.h"
#include "oml/imp/minmax.h"
#include "oml/imp/tstream.h"
#include "oml/imp/cow.h"
#include "oml/random.h"
#include "oml/vector3d.h"
#include <vector>

template <class T> class DiagonalMatrix;


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

enum class FillType {Zero,Random,Unit};

template <class T> class Vector
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
template <class T> Vector<T> operator&(const Vector<T>& a, const Vector<T>& b)
{
  typedef typename Vector<T>::const_iterator CI;
  VecLimits newlim(a.GetLimits().Low,a.GetLimits().Low+a.size()+b.size()-1);
  Vector<T> ret(newlim);
  typename Vector<T>::iterator i=ret.begin();
  for (CI ab=a.begin();ab!=a.end();ab++,i++) *i=*ab;
  for (CI bb=b.begin();bb!=a.end();bb++,i++) *i=*bb;
  return ret;
}

template <class T> void Vector<T>::ReIndex(const std::vector<index_t>& index)
{
  assert(size()==index.size());

  std::vector<index_t>::const_iterator b=index.begin();
  Vector<T>                      dest(GetLimits());
  iterator                       i=dest.begin();
  for (;b!=index.end();b++,i++) *i=(*this)[*b];
  *this=dest;
}

template <class T> Vector<T> Vector<T>::SubVector(const VecLimits& lim) const
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

#undef CHECK

#endif //_Vector_H_
