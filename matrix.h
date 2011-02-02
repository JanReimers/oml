// File: Matrix.H  Matrix class with indirect addressing.
#ifndef _Matrix_H_
#define _Matrix_H_

// Copyright (1994-2005), Jan N. Reimers

#include "oml/iterable.h"
#include "oml/tstream.h"
#include "oml/matrixbase.h"
#include "oml/matsub.h"
#include "oml/matrixalg.h"
#include "oml/cow.h"
#include "oml/matindex.h"
#include "oml/rowcol.h"

template <class T> class Array;

//----------------------------------------------------------------------------
/*! \class Matrix matrix.h oml/matrix.h
  \brief Numerical container with FORTRAN 2-D array symatics.

  Full matrix class with indirect subscripting, which means a list of column pointers is maintained
  to speed up access. This is a standard speed versus excess memory trade off. See DirectMatrix for the other
  option where indexes are calcuated for each access. The importance of the speed difference 
  will depend on your application and CPU/FPU architecture.

  \b Math:
  The special operators for Matrix are
  - \c M2=Transpose(M1) and \c M2=~M1, matrix transpose.
  - \c M \c = \c OuterProduct(V1,V2), std::vector outer product, \c M(i,j)=V1(i)*V2(j) for all \c i,j.
  - \c M3=M1*M2, matrix multiplication.
  - \c V2=M1*V1, matrix-std::vector multiply.
  - \c V2=V1*M1, std::vector-matrix multiply.
  - \c Unit(M), fills \c M(i,j)=delta_ij.
  - \c IsSymmetric(M), returns \c true of \c false.
  - \c IsAntiSymmetric(M), returns \c true of \c false.
  - \c IsHermitian(M), returns \c true of \c false.
  - \c MakeSymmetric(M), returns largest difference between upper and lower elements.
  - \c Normalize(M,V), does \c M(i,j)*=V(i)*V(j) for all \c i,j.
  - \c Normalize(M), does  \c M(i,j)*=1/sqrt(M(i,i)*M(j,j)), for all \c i,j.
  - Eigen routines, matrix inversion, and other numerical operations can be found in \c numeric.h .
  - Complex Eigen routines can be found in \c cnumeric.h

  \b Proxies:
  \c GetRow, \c GetColumn and \c GetDiagonal member functions returns proxies that behave like Vectors. The proxies just
    contain a refernce to the original Matrix. You can even assign to one of the proxies.
    For example:
    \code
    #include "oml/matrix.h"
    Matrix<double> A(10,10),B(10,10);
    // ... file A and B ...
    A.GetRow(2)=sqrt(B.GetDiagonal()); //As usual expression templates generate tight machine code at compile time.
    \endcode

    \nosubgrouping
*/

template <class T> class Matrix 
  : public Indexable<T,Matrix<T>,Full,Real,MatrixShape>
  , public Iterable<T,Matrix<T> >
  , public TStreamableObject<Matrix<T> >
  , public MatrixBase
{
 public:
	
  /*! \name Constructors/Assignment
  Copy constructor and op=(Matrix) are automaically supplied by the compiler.
  */
  //@{
  explicit Matrix(                );//!< Matrix with size=0;
  explicit Matrix(index_t r, index_t c);  //!<Construct from row and column size, use default lower bound.	
  explicit Matrix(subsc_t rl,subsc_t rh,subsc_t cl,subsc_t ch);  //!<Construct from lower and upper bounds.
  explicit Matrix(const VecLimits& r,const VecLimits& c);//!<Construct from lower and upper bounds.
  explicit Matrix(const MatLimits&);//!<Construct from lower and upper bounds.
  //! Copy constructor.
  explicit Matrix(const Matrix& m);
  //! Allows construction from an expression template.
  template <class A,Store M, Data D> Matrix(const Indexable<T,A,M,D,MatrixShape>&);

  //! Assignment.
  Matrix& operator=(const Matrix&);
  //! Allows assignment from an expression template.  
  template <class A,Store M, Data D> Matrix& operator=(const Indexable<T,A,M,D,MatrixShape>&);
  //@}
  
  ~Matrix();

  std::ostream& Write(std::ostream&) const;
  std::istream& Read (std::istream&)      ;
  // Need to disambiguate from expression version.
  friend std::ostream& operator<<(std::ostream& os,const Matrix& a) 
  {
    return os << static_cast<const TStreamableObject<Matrix<T> >& >(a);
  }

  /*! \name Subscripting operators
    FORTRAN style 2-D array subscripting.
    If DEBUG is defined, every index will be checked that it is in range.
    For fast write access with \c op() make a subscriptor, then the COW check is only done once
    during construction of the Subscriptor.
  */
  //@{
  //! const element acces operator, fast and \e cannot trigger a COW operation.
  T  operator()(subsc_t,subsc_t) const;
  //! non-const version can trigger a COW operation, and checks for this with every access.
  T& operator()(subsc_t,subsc_t)      ;
  //@}

  index_t   size     () const; //!<Returns number elements in the Matrix.
  MatLimits GetLimits() const; //!<Returns row/column lower and upper limits a structure.
	
  /*! \name Resize functions
    The data can be optionally preserved.
  */
  //@{
  void SetLimits(const MatLimits&                           , bool preserve=false); //!<Resize from new limits.
  void SetLimits(index_t r, index_t c                       , bool preserve=false); //!<Resize from new rol and column sizes.
  void SetLimits(subsc_t rl,subsc_t rh,subsc_t cl,subsc_t ch, bool preserve=false); //!<Resize from new limits.
  void SetLimits(const VecLimits& r,const VecLimits& c      , bool preserve=false); //!<Resize from new limits.
  //@}

  //! Does \c M(i,j)=M(index[i],j) for i=low...high.  Used for sorting.
  void ReIndexRows   (const Array<index_t>& index) {::ReIndexRows   (*this,index);}
  //! Does \c M(i,j)=M(i,index[j]) for i=low...high.  Used for sorting.
  void ReIndexColumns(const Array<index_t>& index) {::ReIndexColumns(*this,index);}
  //! Does \c M(i,k)=M(j,k) for k=low...high.  Used for sorting.
  void SwapRows   (subsc_t i,subsc_t j) {::SwapRows   (*this,i,j);}
  //! Does \c M(k,i)=M(k,j) for k=low...high.  Used for sorting.
  void SwapColumns(subsc_t i,subsc_t j) {::SwapColumns(*this,i,j);}
  //! Exctract a portion of the matrix. Returns a new Matrix object.
  Matrix SubMatrix(const MatLimits& lim) const {Matrix ret(lim);::SubMatrix(ret,*this);return ret;}

  typedef MatrixRow     <T,Matrix<T>,Full,Real> RowType;
  typedef MatrixColumn  <T,Matrix<T>,Full,Real> ColType;
  typedef MatrixDiagonal<T,Matrix<T>,Full,Real> DiagType;
  
  /*! \name Slice operations.
   These member functions return proxies that behave like Vectors. 
  */
  //@{
  //! Return slice proxy.
  RowType  GetRow     (subsc_t row)       {return RowType (*this,row);} 
  ColType  GetColumn  (subsc_t col)       {return ColType (*this,col);} 
  DiagType GetDiagonal(           )       {return DiagType(*this    );} 
  const RowType  GetRow     (subsc_t row) const {return RowType (*this,row);}
  const ColType  GetColumn  (subsc_t col) const {return ColType (*this,col);}
  const DiagType GetDiagonal(           ) const {return DiagType(*this    );}
  //@}

  /*! \name Iterators.
    Iterators should be STL compatable. These are usefull if you want access all the Matrix
    elements in no particular order.
   */
  //@{
  //! Read only iterator.
  typedef typename Iterable <T,Matrix>::const_iterator  const_iterator ;
  //! Read/write iterator.
  typedef typename Iterable <T,Matrix>::iterator iterator;
  //@}

#if DEBUG
  #define CHECK(i,j) assert(itsLimits.CheckIndex(i,j));
#else
  #define CHECK(i,j)
#endif
//-----------------------------------------------------------------------------
//
//  Allows fast L-value access.  Does COW check during construction.
//
  class Subscriptor 
  {
   public:
    Subscriptor(Indexable<T,Matrix,Full,Real,MatrixShape>& m) 
      : itsLimits(m.GetLimits())
      , itsPtr(static_cast<Matrix&>(m).GetColumns()) 
      {};
		
    T& operator()(subsc_t i,subsc_t j) {CHECK(i,j);return itsPtr[j][i];}
		
   private:
    MatLimits itsLimits;
    T**   itsPtr; 
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
    ArraySubscriptor(Indexable<T,Matrix,Full,Real,MatrixShape>& a) 
      : itsPtr(static_cast<Matrix*>(&a)->Get())
      , itsSize(a.size()) 
      {assert(itsPtr);}
    T& operator[](index_t i) {CHECK(i);return itsPtr[i];}
   private:
    T*      itsPtr;
    index_t itsSize;
  };

#undef CHECK
 
 private:
  friend class Indexable<T,Matrix,Full,Real,MatrixShape>;
  friend class Iterable<T,Matrix>;
  friend class Subscriptor;
  friend class ArraySubscriptor;
	
  T  operator[](index_t i) const;
  T& operator[](index_t i)      ;
//
//  Get pointer to raw data.  required by iterable and Indexable.
//
  const T* Get() const ;       
        T* Get()       ; // Can trigger COW.
  T*const* GetColumns() const {return itsColumns;}
  T**      GetColumns(); // Can trigger COW.
//
//  Reload array of column pointers offset for row lower index.
//
  void NewColumnPointers();  //Reload column pointers.
  void ReleaseColumns();
//
//  Check that column pointers are consistent with the raw data.
//
  bool CheckColumns()       {return &itsColumns[GetLimits().Col.Low][GetLimits().Row.Low]==itsData.Get();} //This should force deep copy.
  bool CheckColumns() const {return &itsColumns[GetLimits().Col.Low][GetLimits().Row.Low]==itsData.Get();} //This should not force deep copy.
//
//  Full internal consistency check.
//
  void          Check () const;       //Check internal consistency between limits and cow.  

  cow_array<T > itsData;     //Copy-On-Write array for the data.
  T**           itsColumns;  //Column pointers.
};

//-----------------------------------------------------------------------------
//
//  Macro, which expands to an index checking function call, 
//  when DEBUG is on
//
#if DEBUG
  #define CHECK(i,j) GetLimits().CheckIndex(i,j)
#else
  #define CHECK(i,j)
#endif

template <class T> inline Matrix<T>::~Matrix()
{
  ReleaseColumns();
}

template <class T> inline T Matrix<T>::operator()(subsc_t i,subsc_t j) const 
{
  CHECK(i,j);
  return GetColumns()[j][i];
}

template <class T> inline T& Matrix<T>::operator()(subsc_t i,subsc_t j)
{
  CHECK(i,j);
  return GetColumns()[j][i]; //Force COW check.
}


#undef CHECK

#if DEBUG
  #define CHECK(i) assert(i>=0 && i<itsData.size())
#else
  #define CHECK(i)
#endif

template <class T> inline  T Matrix<T>::operator[](index_t i) const 
{
  CHECK(i);
  return itsData.Get()[i];
}

template <class T> inline T& Matrix<T>::operator[](index_t i) 
{
  CHECK(i);
  return itsData.Get()[i];
}
#undef CHECK

template <class T> template <class A,Store M, Data D> inline 
Matrix<T>::Matrix(const Indexable<T,A,M,D,MatrixShape>& m)
  : MatrixBase(m.GetLimits        ())
  , itsData   (GetLimits().size())
  , itsColumns(0)
  {
    NewColumnPointers();
    AssignFrom(m); //Use op[i] or op(i,j) depending on M and D.
  }


template <class T> template <class A,Store M, Data D> inline 
Matrix<T>& Matrix<T>::operator=(const Indexable<T,A,M,D,MatrixShape>& m)
{
  if (size()==0) SetLimits(m.GetLimits());
  AssignFrom(m); //Use op[i] or op(i,j) depending on M and D.
  return *this;
}


template <class T> inline index_t Matrix<T>::size() const
{
  return GetLimits().size();
}

template <class T> inline const T* Matrix<T>::Get() const 
{
  return itsData.Get();
}

template <class T> inline void Matrix<T>::SetLimits(const MatLimits& lim,bool preserve)
{
  ::SetLimits(*this,lim,preserve);
}

template <class T> inline void Matrix<T>::SetLimits(index_t r, index_t c , bool preserve)
{
  SetLimits(MatLimits(r,c),preserve);
}

template <class T> inline void Matrix<T>::SetLimits(subsc_t rl,subsc_t rh,subsc_t cl,subsc_t ch, bool preserve)
{
  SetLimits(MatLimits(rl,rh,cl,ch),preserve);
}

template <class T> inline void Matrix<T>::SetLimits(const VecLimits& r,const VecLimits& c , bool preserve)
{
  SetLimits(MatLimits(r,c),preserve);
}

template <class T> inline MatLimits Matrix<T>::GetLimits() const 
{
  return MatrixBase::GetLimits();
}

#endif //_Matrix_H_
