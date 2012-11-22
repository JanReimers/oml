// File: vector2d.h  real space 2 vector.
#ifndef _vector2d_h_
#define _vector2d_h_

#include "oml/indext.h"
#include "oml/mixtypes.h"
#include "oml/Angle.H"

//-----------------------------------------------------------------------------
//
//  Structure for real space two component vectors.  Most algeabraic operators have
//  been overloaded.  So Vector2d's should behave like any intrinsic data type.
//  Some specialized operators are:
//    |  - angle between two vectors in radians.
//    || - angle between two vectors in degrees.
//    !  - magnitude of the vector (assumes an orthogonal isotrpoic metric)
//    ~  - returns unit vector.
//

template <class T> class Vector2d
{
 public:
  Vector2d(                 ): x(  0),y(  0) {};
  Vector2d(T _x,T _y        ): x( _x),y( _y) {};
  Vector2d(const Vector2d& v): x(v.x),y(v.y) {};


  template <class T1> Vector2d(T1 _x,T1 _y          ) : x( _x),y( _y) {}
  template <class T1> Vector2d(const Vector2d<T1>& v) : x(v.x),y(v.y) {}
 ~Vector2d() {};

  Vector2d& operator =(const Vector2d& v) {x=v.x;y=v.y;return *this;}
  template <class T1> Vector2d& operator =(const Vector2d<T1>& v) {x=v.x;y=v.y;return *this;}

  Vector2d<double> GetValue() const {return Vector2d<double>(x.GetValue(),y.GetValue()); }


  const T& operator()(subsc_t i) const {return (&x)[i-1];}
        T& operator()(subsc_t i)       {return (&x)[i-1];}

  T x,y;    //Coordinates.  WOW real data!
};



//-----------------------------------------------------------------------------
//
//  Binary algeabra.
//


template <class T1, class T2> inline
Vector2d<typename BinaryRetType<T1,T2>::RetType> operator +(const Vector2d<T1>& a,const Vector2d<T2>& b)
{
  return Vector2d<typename BinaryRetType<T1,T2>::RetType>(a.x+b.x,a.y+b.y);
}

template <class T1, class T2> inline
Vector2d<typename BinaryRetType<T1,T2>::RetType> operator -(const Vector2d<T1>& a,const Vector2d<T2>& b)
{
  return Vector2d<typename BinaryRetType<T1,T2>::RetType>(a.x-b.x,a.y-b.y);
}

template <class T1, class T2> inline
typename BinaryRetType<T1,T2>::RetType operator *(const Vector2d<T1>& a,const Vector2d<T2>& b)
{
  return (a.x*b.x+a.y*b.y);
}

template <class T> inline
Vector2d<T> operator *(const Vector2d<T>& a,const T F)
{
  return Vector2d<T>(a.x*F,a.y*F);
}

template <class T> inline
Vector2d<T> operator *(const T F,const Vector2d<T>& a)
{
  return Vector2d<T>(a.x*F,a.y*F);
}

template <class T> inline
Vector2d<T> operator /(const Vector2d<T>& a,const T F)
{
  return Vector2d<T>(a.x/F,a.y/F);
}

template <class T1, class T2> inline
Vector2d<typename BinaryRetType<T1,T2>::RetType> operator %(const Vector2d<T1>& a,const Vector2d<T2>& b)
{
  return Vector2d<typename BinaryRetType<T1,T2>::RetType>( a.x%b.x, a.y%b.y );
}

//
//  A operator= B overloads for binary operators
//

template <class T1, class T2> inline
Vector2d<T1>& operator +=(Vector2d<T1>& a,const Vector2d<T2>& b)
{
  a.x+=b.x; a.y+=b.y;
  return a;
}

template <class T1, class T2> inline
Vector2d<T1>& operator -=(Vector2d<T1>& a,const Vector2d<T2>& b)
{
  a.x-=b.x; a.y-=b.y;
  return a;
}

template <class T> inline
Vector2d<T>& operator *=(Vector2d<T>& a,const T F)
{
  a.x*=F;a.y*=F;
  return a;
}

template <class T> inline
Vector2d<T>& operator /=(Vector2d<T>& a,const T F)
{
  a.x/=F;a.y/=F;
  return a;
}

//-----------------------------------------------------------------------------
//
//  Relational operators.
//
template <class T> inline
bool operator ==(const Vector2d<T>& a,const Vector2d<T>& b)
{
  return ((a.x==b.x)&&(a.y==b.y));
}

template <class T> inline
bool operator !=(const Vector2d<T>& a,const Vector2d<T>& b)
{
  return !(a==b);
}

template <class T> inline
bool operator > (const Vector2d<T>& a,const Vector2d<T>& b)
{
  return (!a > !b);
}

template <class T> inline
bool operator < (const Vector2d<T>& a,const Vector2d<T>& b)
{
  return (!a < !b);
}

template <class T> inline
bool operator >=(const Vector2d<T>& a,const Vector2d<T>& b)
{
  return (!a >= !b);
}

template <class T> inline
bool operator <=(const Vector2d<T>& a,const Vector2d<T>& b)
{
  return (!a <= !b);
}

//-----------------------------------------------------------------------------
//
//  Unary operators
//
template <class T> inline
Vector2d<T>  operator -(const Vector2d<T>& a)
{
  return Vector2d<T>(-a.x,-a.y);
}

template <class T> inline
Vector2d<T>  operator +(const Vector2d<T>& a)
{
  return a;
}

//-----------------------------------------------------------------------------
//
//  Angle between two vectors (Radians).
//
template <class T> inline RingAngle operator |(const Vector2d<T>& a,const Vector2d<T>& b)
{
   T phi=fmod(atan2(b.y,b.x)-atan2(a.y,a.x),2*M_PI);
   return phi;
}

template <class T> inline RingAngle arg(const Vector2d<T>& a)
{
   return atan2(a.y,a.x);
}

//-----------------------------------------------------------------------------
//
//  Angle between two vectors (Degrees).
//
/*template <class T> inline
T operator||(const Vector2d<T>& a,const Vector2d<T>& b)
{
  return (T)(acos( (~a) * (~b) )/M_PI*180.0);
  }*/

//-----------------------------------------------------------------------------
//
//  Magnitude and normalize.
//
template <class T> inline
T operator !(const Vector2d<T>& a)
{
  return (T)sqrt(a*a);
}

template <class T> inline
Vector2d<T>  operator ~(const Vector2d<T>& a)  //normalize.
{
  return a/!a;
}

template <class T> inline
Vector2d<T>  FromPolar(T r, RingAngle theta)
{
   return Vector2d<T>(r*cos(theta),r*sin(theta));
}

/*
//--------------------------------------------------------------------
//
//  Specialized templates for complex data types.
//
template <class T> inline Vector2d<complex<T> > conj(const Vector2d<complex<T> >& v)
{
  return Vector2d<complex<T> >(conj(v.x),conj(v.y),conj(v.z));
}

template <class T> inline Vector2d<T> real(const Vector2d<complex<T> >& v)
{
  return Vector2d<T>(real(v.x),real(v.y),real(v.z));
}


template <class T> inline Vector2d<T> imag(const Vector2d<complex<T> >& v)
{
  return Vector2d<T>(imag(v.x),imag(v.y),imag(v.z));
}

template <class T> inline T norm(const Vector2d<complex<T> >& v)
{
  return norm(v.x) + norm(v.y) + norm(v.z);
}
*/


#endif
