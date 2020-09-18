// File: Vector.h  Interface base class for std::vector objects.
#ifndef _Vector_h_
#define _Vector_h_

// Copyright (1994-2005), Jan N. Reimers

#include "oml/array.h"
#include "oml/veclimit.h"
#include "oml/vecindex.h"

template <class T> class complex;
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
template <class T> class Vector
  : public Indexable<T,Vector<T>,Full,Real,VectorShape>
  , public Iterable<T,Vector<T> >
  , public TStreamableObject<Vector<T> >
{
 public:
  /*! \name Constructors/Assignment
  Copy constructor and op=(Vector) are automaically supplied by the compiler.
  */
  //@{
           Vector(                ); //!< Vector with size=0;
  explicit Vector(index_t         ); //!<  All elements are un-initiated, low index is 1.
  explicit Vector(subsc_t,subsc_t ); //!<  Specify lower and upper index.
  Vector(const VecLimits&); //!<  Specify lower and upper index.
  //! Allows construction from an expression template.
  template <class B,Data D> Vector(const Indexable<T,B,Full,D,VectorShape>&);
  //! Allows assignment from an expression template.
  template <class B,Data D> Vector& operator=(const Indexable<T,B,Full,D,VectorShape>&);
  //@}

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
  T  operator()(subsc_t) const;
  //! non-const version can trigger a COW operation, and checks for this with every access.
  T& operator()(subsc_t)      ;
  //@}

  index_t   size     () const; //!<Returns number elements in the Vector.
  VecLimits GetLimits() const; //!<Returns lower and upper limits in a structure.
  subsc_t   GetLow   () const; //!<Returns lower index.
  subsc_t   GetHigh  () const; //!<Returns upper index.

  /*! \name Resize functions
    The data can be optionally preserved.
  */
  //@{
  void SetLimits(const VecLimits&, bool preserve=false); //!<Resize from new limits.
  void SetLimits(index_t         , bool preserve=false); //!<Resize from new size.
  void SetLimits(subsc_t,subsc_t , bool preserve=false); //!<Resize from new limits.
  //@}

  //! Does V(i)=V(index[i]) for i=low...high.  Used for sorting.
  void   ReIndex(const Array<index_t>&);

  /*! \name Sub-vector functions  */
  //@{
  Vector SubVector(const VecLimits&) const; //!< From limits.
  Vector SubVector(index_t         ) const; //!< First n elements.
  Vector SubVector(subsc_t,subsc_t ) const; //!< From limits.
  //@}

  /*! \name Iterators.
    Iterators should be STL compatable.
   */
  //@{
  //! Read only iterator.
  typedef typename Iterable <T,Vector>::const_iterator  const_iterator;
  //! Read/write iterator.
  typedef typename Iterable <T,Vector>::iterator iterator;
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
      , itsPtr(static_cast<Vector*>(&a)->Get()-itsLimits.Low)
      {assert(itsPtr);}

    T& operator()(subsc_t i) {CHECK(i);return itsPtr[i];}

   private:
    VecLimits itsLimits;
    T*        itsPtr;
  };
#undef CHECK

#if DEBUG
  #define CHECK(i) assert(i>=0&&i<itsSize)
#else
  #define CHECK(i)
#endif
  class ArraySubscriptor
  {
   public:
    ArraySubscriptor(Indexable<T,Vector,Full,Real,VectorShape>& a)
      : itsPtr(static_cast<Vector*>(&a)->Get())
      , itsSize(a.size())
      {assert(itsPtr);}
    T& operator[](index_t i) {CHECK(i);return itsPtr[i];}
   private:
    T*      itsPtr;
    index_t itsSize;
  };

#undef CHECK

 private:
  friend class Indexable<T,Vector,Full,Real,VectorShape>;
  friend class Iterable <T,Vector>;
  friend class Subscriptor;
  friend class ArraySubscriptor;
  friend class DiagonalMatrix<T>;

  T  operator[](index_t) const;
  T& operator[](index_t)      ;

  const T* Get() const; //Required by iterable.
        T* Get()      ; //Required by iterable.
  void  Check () const; //Check internal consistency between limits and cow.

  VecLimits    itsLimits; //Manages the upper and lower vector limits.
  cow_array<T> itsData;   //Copy-On-Write array for the data.
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
template <class T, class A, class B, Store MA, Store MB, Data DA, Data DB> inline
T operator*(const Indexable<T,A,MA,DA,VectorShape>& a, const Indexable<T,B,MB,DB,VectorShape>& b)
{
   return Dot(a,b);
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

template <class T> inline T Vector<T>::operator()(subsc_t i) const
{
  CHECK(i);
  return itsData.Get()[itsLimits.Offset(i)];
}

template <class T> inline T& Vector<T>::operator()(subsc_t i)
{
  CHECK(i);
  return itsData.Get()[itsLimits.Offset(i)];
}

#undef CHECK

#if DEBUG
  #define CHECK(i) assert(i>=0 && i<itsData.size())
#else
  #define CHECK(i)
#endif

template <class T> inline  T Vector<T>::operator[](index_t i) const
{
  CHECK(i);
  return itsData.Get()[i];
}

template <class T> inline T& Vector<T>::operator[](index_t i)
{
  CHECK(i);
  return itsData.Get()[i];
}

#undef CHECK

#ifdef DEBUG
#define CHECK Check()
#else
#define CHECK
#endif

template <class T> inline Vector<T>::Vector()
  : itsData(0)
  {
    CHECK;
  }

template <class T> inline Vector<T>::Vector(index_t size)
  : itsLimits(size)
  , itsData  (size)
  {
    CHECK;
  }

template <class T> inline Vector<T>::Vector(index_t l,index_t h)
  : itsLimits(l,h)
  , itsData  (itsLimits.size())
  {
    CHECK;
  }

template <class T> inline Vector<T>::Vector(const VecLimits& lim)
  : itsLimits(lim          )
  , itsData  (lim.size())
  {
    CHECK;
  }

template <class T> inline index_t Vector<T>::size() const
{
  return GetLimits().size();
}

template <class T> inline const T* Vector<T>::Get() const
{
  return itsData.Get();
}

template <class T> inline T* Vector<T>::Get()
{
  return itsData.Get();
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

template <class T> inline subsc_t Vector<T>::GetLow() const
{
  return itsLimits.Low;
}

template <class T> inline subsc_t Vector<T>::GetHigh() const
{
  return itsLimits.High;
}

template <class T> inline void Vector<T>::SetLimits(index_t size,bool preserve)
{
  SetLimits(VecLimits(size),preserve);
}

template <class T> inline void Vector<T>::SetLimits(subsc_t l,subsc_t h,bool preserve)
{
  SetLimits(VecLimits(l,h),preserve);
}

template <class T> inline Vector<T> Vector<T>::SubVector(index_t size) const
{
  return SubVector(VecLimits(size));
}

template <class T> inline Vector<T> Vector<T>::SubVector(subsc_t l,subsc_t h) const
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


#include "oml/minmax.h"

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

      subsc_t low =Max(GetLow() ,theLimits.Low );   //Limits of old and new
      subsc_t high=Min(GetHigh(),theLimits.High);   //data overlap.

      Subscriptor source(*this);         //Destination subscriptor.

      for (subsc_t i=low;i<=high;i++) dest(i)=source(i); //Transfer any overlaping data.
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

template <class T> void Vector<T>::ReIndex(const Array<index_t>& index)
{
  assert(size()==index.size());

  Array<index_t>::const_iterator b=index.begin();
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
  for (subsc_t i=ret.GetLimits().Low;i<=ret.GetLimits().High;i++) r(i)=(*this)(i);
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
