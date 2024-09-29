// File: matrix3d.h
#ifndef _matrix3d_h_
#define _matrix3d_h_

// Copyright (1994-2003), Jan N. Reimers

#include "oml/vector3d.h"

//------------------------------------------------------------------------
//
//  Hard coded 3x3 matrix.  Most algeabraic operators have
//  been overloaded. Also has convenient methods for making
//  3D rotation matricies.
//
template <class T> class Matrix3D
{
public:
   T M11,M12,M13;
   T M21,M22,M23;
   T M31,M32,M33;


  Matrix3D();
  Matrix3D(const T&,const T&,const T&,
 	   const T&,const T&,const T&,
	   const T&,const T&,const T&
	  );
  Matrix3D(T& all);
  Matrix3D(T& diag, T& offdiag);
  Matrix3D(const Matrix3D<T>& m);

  Matrix3D& operator =(const Matrix3D<T>&);
  Vector3D<T> GetRow(index_t) const;
  Vector3D<T> GetCol(index_t) const;

  const T& operator()(index_t i, index_t j) const {return (&M11)[3*(i-1)+j-1];}
        T& operator()(index_t i, index_t j)       {return (&M11)[3*(i-1)+j-1];}
};

//------------------------------------------------------------------------
//
//  Construction zone.
//
template <class T> inline Matrix3D<T>::Matrix3D() :
  M11(1),M12(0),M13(0),
  M21(0),M22(1),M23(0),
  M31(0),M32(0),M33(1)
  {}

template <class T> inline Matrix3D<T>::Matrix3D
  (const T& m11,const T& m12,const T& m13,
   const T& m21,const T& m22,const T& m23,
   const T& m31,const T& m32,const T& m33) :
  M11(m11),M12(m12),M13(m13),
  M21(m21),M22(m22),M23(m23),
  M31(m31),M32(m32),M33(m33)
  {}

template <class T> inline Matrix3D<T>::Matrix3D(T& all) :
  M11(all),M12(all),M13(all),
  M21(all),M22(all),M23(all),
  M31(all),M32(all),M33(all)
  {}

template <class T> inline Matrix3D<T>::Matrix3D(T& diag, T& offdiag) :
  M11(diag   ),M12(offdiag),M13(offdiag),
  M21(offdiag),M22(diag   ),M23(offdiag),
  M31(offdiag),M32(offdiag),M33(diag   )
  {}

template <class T> inline Matrix3D<T>::Matrix3D(const Matrix3D<T>& m) :
  M11(m.M11),M12(m.M12),M13(m.M13),
  M21(m.M21),M22(m.M22),M23(m.M23),
  M31(m.M31),M32(m.M32),M33(m.M33)
  {}

template <class T> inline Matrix3D<T>& Matrix3D<T>::operator=(const Matrix3D<T>& m)
{
  M11=m.M11; M12=m.M12; M13=m.M13;
  M21=m.M21; M22=m.M22; M23=m.M23;
  M31=m.M31; M32=m.M32; M33=m.M33;
  return *this;
}

//------------------------------------------------------------------------
//
//  Comparison
//
template <class T> inline
bool operator==(const Matrix3D<T>& a,const Matrix3D<T>& b)
{
  return a.M11==b.M11 && a.M12==b.M12 && a.M13==b.M13 &&
         a.M21==b.M21 && a.M22==b.M22 && a.M23==b.M23 &&
         a.M31==b.M31 && a.M32==b.M32 && a.M33==b.M33   ;
}
template <class T> inline
bool operator!=(const Matrix3D<T>& a,const Matrix3D<T>& b)
{
	return !(a==b);
}

//------------------------------------------------------------------------
//
//  Container stuff.
//
template <class T> inline Matrix3D<T> Transpose(const Matrix3D<T>& m)
{
	return Matrix3D<T>
     (
      m.M11 , m.M12 , m.M13,
      m.M21 , m.M22 , m.M23,
      m.M31 , m.M32 , m.M33
      );
}

template <class T> inline Matrix3D<T> operator~(const Matrix3D<T>& m)
{
	return Transpose(m);
}

template <class T> inline Vector3D<T> Matrix3D<T>::GetRow(index_t r) const
{
   return Vector3D<T>((*this)(r,1),(*this)(r,2),(*this)(r,3));
}

template <class T> inline Vector3D<T> Matrix3D<T>::GetCol(index_t c) const
{
   return Vector3D<T>((*this)(1,c),(*this)(2,c),(*this)(3,c));
}
//------------------------------------------------------------------------
//
//  Algebra
//
template <class T1, class T2> inline
auto operator+(const Matrix3D<T1>& a,const Matrix3D<T2>& b)
{
  typedef typename ReturnType<T1,T2>::RetType RT;
  return Matrix3D<RT>
  (
    a.M11+b.M11, a.M12+b.M12, a.M13+b.M13,
    a.M21+b.M21, a.M22+b.M22, a.M23+b.M23,
    a.M31+b.M31, a.M32+b.M32, a.M33+b.M33
  );
}

template <class T1, class T2> inline
Matrix3D<T1>& operator+=(Matrix3D<T1>& a,const Matrix3D<T2>& b)
{
  a.M11+=b.M11; a.M12+=b.M12; a.M13+=b.M13;
  a.M21+=b.M21; a.M22+=b.M22; a.M23+=b.M23;
  a.M31+=b.M31; a.M32+=b.M32; a.M33+=b.M33;
  return a;
}

template <class T1, class T2> inline
auto operator-(Matrix3D<T1>& a,const Matrix3D<T2>& b)
{
  typedef typename ReturnType<T1,T2>::RetType RT;
  return Matrix3D<RT>
  (
    a.M11-b.M11, a.M12-b.M12, a.M13-b.M13,
    a.M21-b.M21, a.M22-b.M22, a.M23-b.M23,
    a.M31-b.M31, a.M32-b.M32, a.M33-b.M33
  );
}

template <class T1, class T2> inline
Matrix3D<T1>& operator-=(Matrix3D<T1>& a,const Matrix3D<T2>& b)
{
  a.M11-=b.M11; a.M12-=b.M12; a.M13-=b.M13;
  a.M21-=b.M21; a.M22-=b.M22; a.M23-=b.M23;
  a.M31-=b.M31; a.M32-=b.M32; a.M33-=b.M33;
  return a;
}


template <class T> inline
Matrix3D<T> operator*(const Matrix3D<T>& m, const T& c)
{
   return Matrix3D<T>
     (
      m.M11*c , m.M12*c , m.M13*c,
      m.M21*c , m.M22*c , m.M23*c,
      m.M31*c , m.M32*c , m.M33*c
      );
}

template <class T> inline
Matrix3D<T> operator*(const T& c, const Matrix3D<T>& m)
{
   return m*c;
}

template <class T> inline
Matrix3D<T> operator/(const Matrix3D<T>& m, const T& c)
{
   return Matrix3D<T>
     (
      m.M11/c , m.M12/c , m.M13/c,
      m.M21/c , m.M22/c , m.M23/c,
      m.M31/c , m.M32/c , m.M33/c
      );
}

template <class T> inline T SumSquares(const Matrix3D<T>& a)
{
  return
     a.M11*a.M11 + a.M12*a.M12 + a.M13*a.M13 +
     a.M21*a.M21 + a.M22*a.M22 + a.M23*a.M23 +
     a.M31*a.M31 + a.M32*a.M32 + a.M33*a.M33
     ;
}

template <class T> inline Matrix3D<T> operator-(const Matrix3D<T>& m)
{
   return Matrix3D<T>
     (
      -m.M11 , -m.M21 , -m.M31,
      -m.M12 , -m.M22 , -m.M32,
      -m.M13 , -m.M23 , -m.M33
      );
}

//------------------------------------------------------------------------
//
//  Matrix Algebra
//
template <class T1, class T2> inline
auto operator*(const Matrix3D<T1>& a, const Matrix3D<T2>& b)
{
  return Matrix3D<typename ReturnType<T1,T2>::RetType>
  (
    a.M11*b.M11 + a.M12*b.M21 + a.M13*b.M31,
    a.M11*b.M12 + a.M12*b.M22 + a.M13*b.M32,
    a.M11*b.M13 + a.M12*b.M23 + a.M13*b.M33,
    a.M21*b.M11 + a.M22*b.M21 + a.M23*b.M31,
    a.M21*b.M12 + a.M22*b.M22 + a.M23*b.M32,
    a.M21*b.M13 + a.M22*b.M23 + a.M23*b.M33,
    a.M31*b.M11 + a.M32*b.M21 + a.M33*b.M31,
    a.M31*b.M12 + a.M32*b.M22 + a.M33*b.M32,
    a.M31*b.M13 + a.M32*b.M23 + a.M33*b.M33
  );
}



template <class T1, class T2> inline
auto operator*(const Matrix3D<T2>& m, const Vector3D<T1>& v)
{
   return Vector3D<typename ReturnType<T1,T2>::RetType>
   (
	m.M11*v.x + m.M12*v.y + m.M13*v.z,
	m.M21*v.x + m.M22*v.y + m.M23*v.z,
	m.M31*v.x + m.M32*v.y + m.M33*v.z
   );
}

template <class T1, class T2> inline
auto operator*(const Vector3D<T1>& v, const Matrix3D<T2>& m)
{
   return Vector3D<typename ReturnType<T1,T2>::RetType>
   (
      m.M11*v.x + m.M21*v.y + m.M31*v.z,
      m.M12*v.x + m.M22*v.y + m.M32*v.z,
      m.M13*v.x + m.M23*v.y + m.M33*v.z
   );
}

template <class T> inline T Determinant(const Matrix3D<T>& a)
{
  return
    a.M11*(a.M22*a.M33-a.M32*a.M23) +
    a.M12*(a.M23*a.M31-a.M33*a.M21) +
    a.M13*(a.M21*a.M32-a.M31*a.M22)
  ;
}

template <class T> inline Matrix3D<T> Invert(const Matrix3D<T>& m)
{
  T d=(T)(1.0/Determinant(m));
  return Matrix3D<T>
  (
    (m.M22*m.M33-m.M32*m.M23)*d, (m.M32*m.M13-m.M12*m.M33)*d, (m.M12*m.M23-m.M22*m.M13)*d,
    (m.M31*m.M23-m.M21*m.M33)*d, (m.M11*m.M33-m.M31*m.M13)*d, (m.M21*m.M13-m.M11*m.M23)*d,
    (m.M21*m.M32-m.M31*m.M22)*d, (m.M31*m.M12-m.M11*m.M32)*d, (m.M11*m.M22-m.M21*m.M12)*d
  );
}

//------------------------------------------------------------------------
//
//  Matrix Rotations.
//
template <class T> inline void RotateX(Matrix3D<T>& a,T theta)
{
  a.M12=a.M13=a.M21=a.M31=0;
  a.M11=1;
  a.M32=-(a.M23=(T)sin(theta));
  a.M22= a.M33=(T)cos(theta);
}

template <class T> inline void RotationY(Matrix3D<T>& a,T theta)
{
  a.M21=a.M23=a.M12=a.M32=0;
  a.M22=1;
  a.M31=-(a.M13=(T)sin(theta));
  a.M11= a.M33=(T)cos(theta);
}

template <class T> inline void RotationZ(Matrix3D<T>& a,T theta)
{
  a.M13=a.M23=a.M31=a.M32=0;
  a.M33=1;
  a.M21=-(a.M12=(T)sin(theta));
  a.M11= a.M22=(T)cos(theta);
}



#endif //_Matrix3D_h_
