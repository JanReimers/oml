// File: vector3d.h  real space 3 vector.
#ifndef _vector3d_h_
#define _vector3d_h_

// Copyright (1994-2003), Jan N. Reimers

#include "oml/indext.h"
#include "oml/mixtypes.h"
#include <cmath>
namespace std {template <class T> class complex;}

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
template <class T> class Vector3D
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
  const T& operator()(subsc_t i) const {return (&x)[i-1];}
  //! Element access
  T& operator()(subsc_t i)       {return (&x)[i-1];}
   

  /*! \name Coordinates*/
  //@{
  T x; //!< \a x coordinate.
  T y; //!< \a y coordinate.
  T z; //!< \a z coordinate. 
  //@}
};


//-----------------------------------------------------------------------------
//
//  Binary algeabra.
//


template <class T1, class T2> inline
Vector3D<typename BinaryRetType<T1,T2>::RetType> operator +(const Vector3D<T1>& a,const Vector3D<T2>& b)             
{
  return Vector3D<typename BinaryRetType<T1,T2>::RetType>(a.x+b.x,a.y+b.y,a.z+b.z);
}

template <class T1, class T2> inline
Vector3D<typename BinaryRetType<T1,T2>::RetType> operator -(const Vector3D<T1>& a,const Vector3D<T2>& b)
{
  return Vector3D<typename BinaryRetType<T1,T2>::RetType>(a.x-b.x,a.y-b.y,a.z-b.z);
}

template <class T1, class T2> inline 
typename BinaryRetType<T1,T2>::RetType operator *(const Vector3D<T1>& a,const Vector3D<T2>& b)
{
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
Vector3D<typename BinaryRetType<T1,T2>::RetType> operator %(const Vector3D<T1>& a,const Vector3D<T2>& b)
{
  return Vector3D<typename BinaryRetType<T1,T2>::RetType>( a.x%b.x, a.y%b.y, a.z%b.z );
}

template <class T1, class T2> inline
Vector3D<typename BinaryRetType<T1,T2>::RetType> Cross(const Vector3D<T1>& a,const Vector3D<T2>& b)
{
  return Vector3D<typename BinaryRetType<T1,T2>::RetType>( a.y*b.z-a.z*b.y, a.z*b.x-a.x*b.z, a.x*b.y-a.y*b.x );
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
  return (!a > !b);
}

template <class T> inline
bool operator < (const Vector3D<T>& a,const Vector3D<T>& b)
{
  return (!a < !b);
}

template <class T> inline
bool operator >=(const Vector3D<T>& a,const Vector3D<T>& b)
{
  return (!a >= !b);
}

template <class T> inline
bool operator <=(const Vector3D<T>& a,const Vector3D<T>& b)
{
  return (!a <= !b);
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
T operator |(const Vector3D<T>& a,const Vector3D<T>& b)
{
  return (T)acos( (~a) * (~b) );
}

//-----------------------------------------------------------------------------
//
//  Angle between two vectors (Degrees).
//
template <class T> inline 
T operator||(const Vector3D<T>& a,const Vector3D<T>& b)
{
  return static_cast<T>(acos( (~a) * (~b) )/M_PI*180.0);
}

//-----------------------------------------------------------------------------
//
//  Magnitude and normalize.
//
template <class T> inline 
T operator !(const Vector3D<T>& a)
{
  return (T)sqrt(a*a);
}

template <class T> inline 
Vector3D<T>  operator ~(const Vector3D<T>& a)  //normalize.
{
  return a/!a;
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

template <class T> inline T norm(const Vector3D<std::complex<T> >& v)
{
  return norm(v.x) + norm(v.y) + norm(v.z);
}


#endif
