#ifndef FAKEDOUBLE_H_INCLUDED
#define FAKEDOUBLE_H_INCLUDED

//
//  Make some fake, no-op for double and double containers, so we can use the same template code
//  for complex and double types
//
inline const double& conj(const double& d) { return d;}
inline const double& real(const double& d) { return d;}
inline       double  imag(const double& d) { return 0.0;}
inline       double  Max (const double& d) { return d;}
//using std::real;

template <class T> class Matrix;
template <class T> class Vector;
template <class T> class Vector3D;

inline const Matrix  <double>& conj(const Matrix  <double>& m) {return m;}
inline const Vector  <double>& conj(const Vector  <double>& v) {return v;}
inline const Vector3D<double>& conj(const Vector3D<double>& v) {return v;}


#endif // FAKEDOUBLE_H_INCLUDED
