// File: Array.h  Array class for any numerical data type.
#ifndef _Array_h_
#define _Array_h_

// Copyright (1994-2005), Jan N. Reimers

#include "oml/indext.h"
#include "oml/shape.h"
#include "oml/tstream.h"
#include "oml/iterable.h"
#include "oml/cow.h"
#include "oml/arrindex.h"

/*! \class Array array.h oml/array.h
  \brief Numerical container with C array symmantics.

  Arrays have element indexes ranging from 0...n-1 and are indexed using the \c [i]
  syntax, just like C arrays.  The base index 0 cannot be changed. Arrays support
  copy on write (COW) symmantics, which avoids deep copying (element by element)
  when they are returned by value from a function. 

  \b Math:

  The special operators for Array are.
  - \c A*A does a direct multiply (not a dot product).
  - use \c Dot(A,A) to get a dot product.
  - \c A/A does a direct divide.

  \nosubgrouping
*/
template <class T> class Array
  : public Indexable<T,Array<T>,Full,Real,ArrayShape>
  , public Iterable<T,Array<T> >
  , public TStreamableObject<Array<T> >
{
 public:

  /*! \name Constructors/Assignment
  Copy constructor and op=(Array) are automatically supplied by the compiler.

  */
  //@{
           Array(              ) : itsData(0        ) {}; //!< Array with size=0;
  explicit Array(index_t   size) : itsData(size     ) {}; //!< All elements are un-initiated.
  //! Allows construction from an expression template.
  template <class B> Array(const Indexable<T,B,Full,Real,ArrayShape>&);
  //! Allows assignment from an expression template.
  template <class B> Array& operator =(const Indexable<T,B,Full,Real,ArrayShape>&);
  //@}


  std::ostream& Write(std::ostream&) const; 
  std::istream& Read (std::istream&)      ; 
  // Need to disambiguate from expression version.
  friend std::ostream& operator<<(std::ostream& os,const Array& a) 
  {
    return os << static_cast<const TStreamableObject<Array<T> >& >(a);
  }

  /*! \name Subscripting operators
    If \c DEBUG is defined, every index will be checked that it is in range.
    For fast write access with \c op[] make a \c Subscriptor, then the COW check is only done once
    during construction of the Subscriptor.
  */
  //@{
  //! const element access operator, fast and \e cannot trigger a COW operation.
  T  operator[](index_t i) const;
  //! non-const version can trigger a COW operation, and checks for this with every access.
  T& operator[](index_t i)      ;
  //@}

  //! Returns number elements in the array. Should be the same as the allocated space.
  index_t size() const;
  //! Change the size of the array and optionally preserve as many values as possible.
  void    SetSize(index_t, bool preserve=false)      ;

  //! Does A[i]=A[index[i]] for i=0...n-1.  Used for sorting.
  void  ReIndex(const Array<index_t>& index);
  //! Extract a subarray.  Involves a deep copy.
  Array SubArray(index_t start, index_t stop) const;

  
#if DEBUG
  #define CHECK(i) assert(i>=0&&i<itsSize)
#else
  #define CHECK(i)
#endif
  class ArraySubscriptor 
  {
   public:
    ArraySubscriptor(Indexable<T,Array,Full,Real,ArrayShape>& a) 
      : itsPtr(static_cast<Array&>(a).Get()), itsSize(a.size()) {assert(itsPtr);}
    T& operator[](index_t i) {CHECK(i);return itsPtr[i];}
   private:
    T*      itsPtr;
    index_t itsSize;
  };
#undef CHECK
  
  /*! \name Subscriptors and Iterators.
    Iterators should be STL compatible.
   */
  //@{
  //! Read only iterator.
  typedef typename Iterable <T,Array>::const_iterator  const_iterator ;
  //! Read/write iterator.
  typedef typename Iterable <T,Array>::iterator iterator;
  //! Default Subscriptor.
  typedef ArraySubscriptor Subscriptor;
  //@}
 private:
  friend class ArraySubscriptor;
  friend class Iterable <T,Array>;
  
  const T* Get() const;
        T* Get()      ;
  void  CheckIndex(index_t) const;

  cow_array<T> itsData;
};


//-------------------------------------------------------------------------------
//
//  These macros invoke Array index bounds checking if DEBUG is on.
//
#if DEBUG
  #define CHECK(i) CheckIndex(i)
#else
  #define CHECK(i)
#endif

template <class T> inline  T Array<T>::operator[](index_t i) const 
{
  CHECK(i);
  return itsData.Get()[i];
}

template <class T> inline T& Array<T>::operator[](index_t i) 
{
  CHECK(i);
  return itsData.Get()[i];
}

#undef CHECK

//-------------------------------------------------------------------------------
//
//  Some inline functions.
//
template <class T> inline index_t  Array<T>::size() const
{ 
  return itsData.size();
}

template <class T> inline const T* Array<T>::Get() const 
{
  return itsData.Get();
}

template <class T> inline T* Array<T>::Get()
{
  return itsData.Get();
}

template <class T> template <class B> inline 
Array<T>& Array<T>::operator=(const Indexable<T,B,Full,Real,ArrayShape>& a) 
{
  ArrayAssign(*this,a);
  return *this;
}

template <class T> template <class B> inline 
Array<T>::Array(const Indexable<T,B,Full,Real,ArrayShape>& a) 
  : itsData(a.size())
  {
    ArrayAssign(*this,a);
  }

//------------------------------------------------------------------------------
//
//  Change the array to a new size.  As much old data as possible
//  is saved if the \param preserve flag is set.
//
template <class T> inline void Array<T>::SetSize(index_t newsize, bool preserve)
{
  assert(newsize >= 0);
  if (size()!=newsize)
  {
    if (preserve)
    {
      Array<T> dest(newsize);
      const_iterator  b=this->begin();
      iterator        i=dest.begin();
      for (;b!=this->end()&&i!=dest.end();b++,i++) *i=*b;
      *this=dest;
    }
     else
    {
      *this=Array<T>(newsize);
    }
  }
}

#if DEBUG
//------------------------------------------------------------------------------
//
//  Check if an index is within bounds.
//
template <class T> inline void Array<T>::CheckIndex(index_t i) const
{
  if (i>=size() || i<0)
  {
    OMLArrayIndexError(i,size()); //Call non inline routine to print the message.
    assert(i<size()); //Leave asserts here so users sees context.
    assert(i>=0);
  }
}
#endif

//! Concatenation of two arrays.
template <class T> inline Array<T> operator&(const Array<T>& a, const Array<T>& b)
{
  typedef typename Array<T>::iterator I;
  typedef typename Array<T>::const_iterator CI;
  Array<T> ret(a.size()+b.size());
  I i=ret.begin();
  for (CI ab=a.begin();ab!=a.end();ab++,i++) *i=*ab;
  for (CI bb=b.begin();bb!=b.end();bb++,i++) *i=*bb;
  return ret;
}

template <class T> inline void Array<T>::ReIndex(const Array<index_t>& index)
{
  assert(size()==index.size());
  
  Array<index_t>::const_iterator b=index.begin();
  Array<T>                       dest(size());
  iterator                       i=dest.begin();
  for (;b!=index.end();b++,i++) *i=(*this)[*b];
  *this=dest;
}

template <class T> inline Array<T> Array<T>::SubArray(index_t start, index_t stop) const
{
  assert(start>=0);
  assert(stop>=start);
  assert(stop<size());
  Array<T> ret(stop-start+1);
  iterator r=ret.begin();
  for (index_t i=start;r!=ret.end();i++,r++) *r=(*this)[i];
  return ret;
}

/*! \mainpage Object %Matrix Library (OML)
 
  \section motive Motivation
  One of the strong points of the FORTRAN language for numerical programming is its support
  for multidimensional arrays. Beyond this FORTRAN (at least the 77 standard) has a 
  number of problems:
  -# No support for dynamic memory allocation. The maximum problem size has to decided at 
  compile time.
  -# No concept of user defined data types and arrays thereof.
  -# Very poor static (compile time) checking of function interfaces etc.
  -# No support for separation of interface and implementation, i.e. no way to do 
  dependency inversions (ref. Martin).
  -# Very poor support for std::string handling and system level operations.

  The C language solves some of these problems but introduces more:
  -# Arrays start from 0 instead of 1. 
  -# Dynamically allocated multi dimensional arrays are hideously complex and error prone.
  -# No integer power operator like A**3.
  -# Static arrays are stored in row major form, this is ass backwards for FORTRAN programmers.

  The bare C++ language again solves some of the above problems but not all. However C++ is
  an extremely expressive language and just about anything can be added to the language
  by means of user defined types and templates.

  \section oml OML Introduction

  First some nomenclature: Instead of the term \e array I will use the term \e container to 
  refer to any C++ entity that can store and manage multiple values of a single type.
  The OML attempts to address most of the FORTRAN and C issues listed above, by providing the following features:
  - 1-D and 2-D containers that are dynamically allocated behind the scenes. So users rarely
  have to think about memory management issues.
  - Containers can be freely passed by value (rather than by reference) with very little 
  speed penalty. This is achieved by simply transferring a data pointer into the
  receiving container. This is called a \e shallow \e copy. The problem of how to deal with a 
  subsequent write operation on a container sharing data with others is addressed below.
  - Containers with user selectable lower indexes that default 1 are defined. This emulates
  FORTRAN.
  - Math and matrix/std::vector algebra functions and operators are defined for all containers. This
  removes the need for many of the "do loops" that tend to propagate throughout numerical
  programs. Example 1:
  \code
  // Instead of writing
  for (int i=1;i<=n;i++) A(i)=sin(B(i));
  // You just write
  A=sin(B);
  \endcode
  Example 2:
  \code
  // Instead of writing
  for (int i=1;i<=n;i++) 
    for (int j=1;j<=n;j++)
    {
      C(i,j)=0.0;
      for (int k=1;k<=n;k++)
        C(i,j)+=A(i,k)*B(k,j);
     }
  // You just write
  C=A*B; 
  \endcode
  
  \section containers OML Containers
  All OML containers are templated so you can put any (hopefully numerical) data type
  in the containers.  OML provides the following containers:
  - \c Array<T> which emulates a 1-D C style array.  Lower index is always 0, and elements 
  are accessed using \c A[i] notation.
  - \c List<T> is like %Array<T> but will dynamically resize itself when the user pushes
  values onto the back of the list (just like and STL \c std::vector<T>).
  - \c Vector<T> which emulates a 1-D FORTRAN array. Lower index is selectable, defaulting 
  to 1. Elements are accessed using \c V(i) notation.
  - \c Matrix<T> which emulates a 2-D FORTRAN array. Lower indices are selectable, defaulting 
  to 1. Elements are accessed using \c M(i,j) notation.
  - \c SMatrix<T> which emulates a 2-D FORTRAN array but only the upper triangle is stored.
  So assigning to \c S(i,j) implicitly assigns \c S(j,i), hence \c S(i,j)==S(j,i) 
  is \b always true. This type of matrix saves some memory but you pay for this in performance
  when accessing the elements.
  - \c SMatrix<std::complex<T> > will magically behave like a Hermitian matrix, so that \c S(i,j)==conj(S(j,i))
  is always true.

  That's it! There are no 3-D or 7-D containers as in FORTRAN. I find that with C++ coding 
  styles you don't need these higher dimensional containers. You could declare a 
  \c Vector<Matrix<double> > or a \c Matrix<Vector<double> > if you are stuck.

  \section type Recommended Data Types.
  So what types can you put in these containers? In principal anything can go in there
  with the following caveats:
  -# If you try something non-sensicle:
  \code
  #include "oml/array,h"
  #include "Froggy.H" //Some user defined class.
  Array<Froggy> froggies(1000); //Ok so far.
  Array<double> B=sin(froggies); //Compile time error, no function sin(Froggy) defined.
  \endcode
  If the user did happen to define a function called \c sin() that takes a Froggy argument
  and returns a double, then the example above will indeed compile. 
  -# The \c Array<T>, \c List<T> and \c Vector<T> are all 100% defined in the header files
  so the user can instance these with any data type. The 2-D matrix containers have a
  fair amount of their code not in the headers. This \e non-inline code is compiled into
  the library liboml.a (or liboml.so) for types T=double and std::complex<double>. So filling
  matrices with other data types will require compiling extra files to avoid linker errors.
  -# Some of the member functions for \c Array<T>, \c List<T> and \c Vector<T> are not
  inline. These are functions that tend to be rarely used such as \c SetSize(). Nevertheless
  there may be some code bloat observed if these containers are instanced for a large
  number of different data types. In most cases it makes more sense to use STL containers
  for storing (by value) non numerical data types.

  In summary I recommend storing only numerical data types in the OML containers. This
  will allow one to take advantage of all overloaded math functions and operators.

  \section memory Memory and Speed Details

  All OML containers support delayed
  copy on write (COW) symmantics, which avoids deep copying (element by element)
  when they are returned by value from a function. Example:
  \code
  #include "oml/array.h"
  Array<double> A(100000);
  // ... fill A ...
  Array<double> B(A); //Does shallow copy, i.e. just transfers the double* from A to B
  //Now B and A share the same data.
  B[0]=1.0; //Any write operation to either A or B triggers a deep Copy On Write(COW).
  //B and A now have their own data arrays.  100000 elements were copied, rather expensive!
  A[0]=2.0; //No COW now because A no longer shares data with B (or any other Array).
  \endcode

  This COW system is also often referred to \e Lazy \e Copying or \e Delayed \e Copying.
  The main reason to implement this is to allow functions to construct containers and return
  them by value without a speed penalty. Example:
  \code
  Array<double> Frobnicate(const Array<double>& A,const Array<double>& B)
  {
    Array<double> C=(A-B)*(A+B)+B*log(A)-exp(A*B)*A+cos(sin(tan(A/b)));
    return C; //this will return a shallow copy of C so very little speed penalty.
  }
  \endcode

  If you return a pointer to the new container 
  the user has to worry about memory management which is a no no. Also returning a reference
  to a temporary (such as \c C) is illegal.

  This \e shallow \e copy system has a penalty. Any access to the individual elements of
  a non constant OML container must first check if there are multiple owners.  I call 
  this a COW check, and it can give a performance hit. This can be solved by using a 
  \c Subscriptor helper class:
  \code
  #include "oml/array.h"
  #include <cmath>

  const double Pi=acos(-1.0);
  int N=1000;

  Array<double> A(N);
  for (int i=0;i<N;i++) A[i]=sin(2*Pi*i*0.01); // Inefficient, does COW check for every access.

  Array<double>::Subscriptor s(A); //Only one COW check is done here.
  for (int i=0;i<N;i++) s[i]=sin(2*Pi*i*0.01); // Faster, but you may have to measure 
                                               // carefully to notice.
  \endcode
  Beware that even reading from \c A in the above example will also invoke a COW check for every access.
  Why? Because \c A is not a constant object. The way to solve this is to froce  the compiler to use
  the const version of \c op[]:
  \code
  #include "oml/array.h"
  Array<double> A(10);
  // ... fill A ...
  const Array<double>& r(A); //No COW check, don't need it.
  double a0=r[0]; //Reference is constant, calls op[] const, so no COW check or deep copy.
  \endcode

  Beware of the following:
  \code
  #include "oml/array.h"
  Array<double> A(10),B(10);
  // ... fill A and B ...
  Array<double>::Subscriptor s(A); //Only one COW check is done here.
  B=A; //Shallow copy, so B and A share the same data.
  s[0]=1.0; //woops!
  \endcode
  We faked out the system, \c 1.0 get assigned to \c A[0] \e and to \c B[0] !!!!!


  Now one can ask about the speed penalty associated with nesting function calls, for example:
  \code
  #include "oml/array.h"
  Array<double> A(10000),B(00010),C(10000);
  // ... fill A and B ...
  C=Dot(Add(A,B),Subtract(A,B)); //Trying to calculate (A+B) dot (A-B).
  \endcode
  There will be a sequence of intermediate %Array's created in this example, and 3 separate
  loops will be executed. This can also be solved. Expression templates are used to convert arbitrary
  algebraic expressions into tight loops, with no creation of intervening temporary
  objects. Example:
  \code
  #include "oml/array.h"
  int N=1000;
  Array<double> A(N),B(N);
  // ... fill A and B ...
  Array<double> C=(A-B)*(A-B)+0.5*(A+B); //Expression template created to build C from.
  \endcode

  The algebraic expression above gets translated at compile time into the following:
  \code
  for (int i=0;i<N;i++) C[i]=(A[i]-B[i])*(A[i]-B[i])+0.5*(A[i]+B[i]);
  \endcode
  The optimizer will then do common sub-expression elimination etc.  Basically you
  get full speed at run time and no creation of extra temporaries. 

  One problem with 
  expression templates is that if there is a compilation error, the error messages are very hard
  read because of the many levels of template nesting.

  \section io IO and Streaming
  OML containers are streamable, so you can write them to a file in ascii or 
  binary form using \c op<<, and then read them back in using \c op>>. The object type
  can even be automatically detected and created.  This is sometimes referred to as 
  \e pickling and \e unpickling!
  Output of an OML object consists of
  -# flag indicating the mode (binary,ascii or pretty)
  -# typeinfo for the object
  -# size and/or limits of the container
  -# all the data in one unindexed list
  
  When OML objects are read in using op>> the following occurs
  -# The mode is set to mode used to save the object.
  -# The typeinfo is checked for a match.
  -# The size/limits info is read in and the object is resized accordingly. Old data is \b not preserved.
  -# The data is read in a stored in the object.
  -# The mode is restored to its original state.
  
  Example:
  \code
  #include "oml/array.h"
  #include "oml/array_io.h"
  #include <iostream>
  #include <fstream>

  StreamableObject::SetToBinary(); //Most compact form.
  Array<double> A(100000);
  // ... fill A ...
  std::cout << A << std::endl;
  ofstream out("test.dat"):
  out << A;
  ifstream in("test.dat");
  Array<double> B; //Don't need to know the size.
  in >> B; 
  if (A==B)
    std::cout << "OK" << std::endl;
  else
    std::cout << "Something wrong!!!" << std::endl;
  \endcode

  \section math Math Function and Operators
  OML defines a large number for math functions and overloaded operators that
  are common to all OML containers.  Everything is inline, and functions can be
  nested.

  The following operators are supported for all containers:
  \code
  +,-,*,/,%, unary + and unary -, ==, !=
  \endcode
  Operator \c * has a different meaning depending the container type. These operators can
  be mixed with the scalers in most cases.
  
  The following math functions are also supported for all containers:
  \code
  DirectMultiply, DirectDivide, Dot, Sum, Min, Max, Integrate
  pow, fmod, atan2, hypot, Equal
  sin, cos, tan, asin, acos, atan, sinh, cosh, tanh, 
  exp, log, pow2, pow3, pow4, pow10, log10, sqrt, fabs, conj, abs, arg, norm, real, imag.
  \endcode
  Some functions are specific to containers of complex values, so misuse will result in
  a compile time error.
   
  Note that multiplication of
  %Array's does not imply an algebraic dot product, use \c Dot(A,B) for that.

  \todo Construct from \c T[] to allow for example \c Array<double> \c A={1.,2.,3.,4.};  
  \todo Construct and assign from STL iterators?

  \author Jan Reimers
  \date Original version Dec. 1992
  \date Put in expression templates Sept. 1995
  \date Major upgrades Feb. 2002
 */



#endif //_Array_h_
