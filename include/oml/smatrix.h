// File: SMatrix.H  Matrix class with symmteric storage.
#ifndef _SMatrix_H_
#define _SMatrix_H_

// Copyright (1994-2005), Jan N. Reimers

#include "oml/matrix.h" //Need some operator definitions.
#include "oml/imp/matrixbase.h"
#include "oml/imp/arrindex.h"
#include "oml/imp/smatindex.h"
#include "oml/imp/tstream.h"
#include "oml/imp/rowcol.h"
#include "oml/imp/cow.h"
#include <vector>

//----------------------------------------------------------------------------
/*! \class SMatrix smatrix.h oml/smatrix.h
  \brief Numerical container symmetric matrix storage symmantics.

  Symmetric matrix class with direct subscripting, which means not list of column pointers is maintained.
  Offsets into the data must calculated for each access, rather expensive in terms of CPU time.
  For this container, saving memory takes priority over speed.
  Most operations supported by Matrix are also supported here.
  An SMatrix<std::complex<double> > will automatically be Hermitian, i.e. S(i,j)=conj(S(j,i))!  This is done
  through the majic of template spcialization.

  \b Math:
  The special operators for Matrix are
  - \c S \c = \c OuterProduct(V), std::vector outer product, \c S(i,j)=V(i)*V(j) for all \c i,j.
  - \c S3=S1*S2, \c M2=S*M1, \c M2=M1*S, matrix multiplication.
  - \c V2=S*V1, matrix-std::vector multiply.
  - \c V2=V1*S, std::vector-matrix multiply.
  - \c Unit(S), fills \c S(i,j)=delta_ij.
  - \c Normalize(S,V), does \c S(i,j)*=V(i)*V(j) for all \c i,j.
  - \c Normalize(S), does  \c S(i,j)*=1/sqrt(S(i,i)*S(j,j)), for all \c i,j.
  - Eigen routines, matrix inversion, and other numerical operations can be found in \c numeric.h .
  - Complex Eigen routines can be found in \c cnumeric.h

  \nosubgrouping
*/
template <class T> class SMatrix
  : public MatrixBase
  , public ArrayIndexable<T,SMatrix<T>,Symmetric     ,MatrixShape>
  , public      Indexable<T,SMatrix<T>,Symmetric,Real,MatrixShape>
  , public TStreamableObject<SMatrix<T> >
{
 public:
  typedef ArrayIndexable<T,SMatrix<T>,Symmetric     ,MatrixShape> ArrayIndexT;
  typedef      Indexable<T,SMatrix<T>,Symmetric,Real,MatrixShape> IndexableT;
  typedef Ref<T,IndexableT,MatrixShape> RefT;

  /*! \name Constructors/Assignment
  Copy constructor and op=(Vector) are automaically supplied by the compiler.
  */
  //@{
  explicit SMatrix(                    ); //!< Matrix with size=0;
  explicit SMatrix(index_t); //!< Construct from row and same column size, use default lower bound.
  explicit SMatrix(const VecLimits& r,const VecLimits& c); //!< For compatability with Matrix.
  explicit SMatrix(const MatLimits&); //!< For compatability with Matrix.
  //! Copy constructor.
  explicit SMatrix(const SMatrix& m);
  //! Allows construction from an approximately sym matrix stored as a full matrix.
  SMatrix(const Matrix<T>&,double tol=0.0);
  //! Allows construction from an expression template.
  template <class B, Data D> SMatrix(const Indexable<T,B,Symmetric,D,MatrixShape>&);

  //! Assignment.
  SMatrix& operator=(const SMatrix&);
  //! Allows assignment from an expression template.
  template <class B, Data D> SMatrix& operator= (const Indexable<T,B,Symmetric,D,MatrixShape>&);
  //@}

//  using Indexable<T,SMatrix<T>,Symmetric,Real,MatrixShape>::AssignFrom;

  std::ostream& Write(std::ostream&) const;
  std::istream& Read (std::istream&)      ;
  // Need to disambiguate from expression version.
  friend std::ostream& operator<<(std::ostream& os,const SMatrix& a)
  {
    return os << static_cast<const TStreamableObject<SMatrix<T> >& >(a);
  }

  /*! \name Subscripting operators
    FORTRAN style 2-D array subscripting. These are slow!
    If DEBUG is defined, every index will be checked that it is in range.
    For fast write access with op() make a subscriptor, then the COW check is only done once
    during construction of the Subscriptor.
  */
  //@{
  //! const element acces operator, fast and \e cannot trigger a COW operation.
  //T  operator()(index_t i,index_t j) const {return SymmetricSubscriptor<T>::const_apply(i,j,GetLimits(),priv_begin());}
  //! non-const version can trigger a COW operation, and checks for this with every access.
  const T& ref(index_t i,index_t j) const;
  T operator()(index_t i,index_t j) const;
//  {
//      index_t index=GetSymOffset(i,j,GetLimits());
//      return i<=j ? itsData[index] : conj(itsData[index]);
//  }
  T& operator()(index_t i,index_t j);
//  {
//      assert(i<=j);
//      index_t index=GetSymOffset(i,j,GetLimits());
//      return itsData[index];
//  }
  //@}

  size_t    size     () const; //!<Returns number elements in the Matrix.
  MatLimits GetLimits() const; //!<Returns row/column lower and upper limits a structure.

  /*! \name Resize functions
    The data can be optionally preserved.
    Many of these require redundant information for compatability with Matrix.
  */
  //@{
  //! Resize and optionally preserve as much data as possible.
  void SetLimits(size_t N       , bool preserve=false);
  void SetLimits(const VecLimits&, bool preserve=false);
  void SetLimits(const MatLimits&, bool preserve=false);
  //@}

  //! Does M(i,j)=M(index[i],j) for i=low...high.  Used for sorting.
  void ReIndexRows   (const std::vector<index_t>& index);
  //! Does M(i,j)=M(i,index[j]) for i=low...high.  Used for sorting.
  void ReIndexColumns(const std::vector<index_t>& index);
  //! Does M(i,k)=M(j,k) for k=low...high.  Used for sorting.
  void SwapRows   (index_t i,index_t j);
  //! Does M(k,i)=M(k,j) for k=low...high.  Used for sorting.
  void SwapColumns(index_t i,index_t j);
  //! Exctract a portion of the matrix. Returns a new Matrix object.
  SMatrix SubMatrix(const MatLimits& lim) const;

  typedef MatrixRow     <T,SMatrix,Symmetric,Real> RowType;
  typedef MatrixColumn  <T,SMatrix,Symmetric,Real> ColType;
  typedef MatrixDiagonal<T,SMatrix,Symmetric,Real> DiagType;

  /*! \name Slice operations.
   These member functions return proxies that behave like Vectors.
  */
  //@{
  //! Return slice proxy.
        RowType  GetRow     (index_t row)       {return RowType (*this,row);}
        ColType  GetColumn  (index_t col)       {return ColType (*this,col);}
        DiagType GetDiagonal(           )       {return DiagType(*this    );}
  const RowType  GetRow     (index_t row) const {return RowType (*this,row);}
  const ColType  GetColumn  (index_t col) const {return ColType (*this,col);}
  const DiagType GetDiagonal(           ) const {return DiagType(*this    );}
  //@}

  /*! \name Iterators.
    Iterators should be STL compatable. These are usefull if you want access all the Matrix
    elements in no particular order.
   */
  //@{
  //! Read only iterator.
  typedef typename ArrayIndexT::const_iterator  const_iterator ;
  //! Read/write iterator.
  typedef typename ArrayIndexT::iterator iterator;
  //@}

  static index_t GetSymOffset(index_t i, index_t j, const MatLimits& lim)
  {
      return i<j ? GetSymOffset(i-lim.Row.Low,j-lim.Col.Low,lim.GetNumRows())
      : GetSymOffset(j-lim.Row.Low,i-lim.Col.Low,lim.GetNumRows());
  }

  static index_t GetSymOffset(index_t i, index_t j, index_t n)
  {
#ifdef UPPER_ONLY
		assert(j>=i);
#endif
#ifdef LOWER_ONLY
		assert(j<=i);
#endif
    return i<j ? (i*(2*n-i-1))/2+j : (j*(2*n-j-1))/2+i;
  }

#if DEBUG
  #define CHECK(i,j) assert(itsLimits.CheckIndex(i,j));
#else
  #define CHECK(i,j)
#endif

  SMatrix& operator+=(const ArrayIndexable<T,SMatrix,Symmetric,MatrixShape>& b) {return ArrayAdd(*this,b);}
  template <class B> SMatrix& operator+=(const Indexable<T,B,Symmetric,Abstract,MatrixShape>& b)
  {
      if (size()==0)
      {
        SetLimits(b.GetLimits(),false);
        Fill(*this,T(0));
      }
      return MatrixAdd(*this,b);
  }
//-----------------------------------------------------------------------------
//
//  Allows fast L-value access.  Does COW check during construction.
//
  class Subscriptor
  {
  public:
    Subscriptor(Indexable<T,SMatrix,Symmetric,Real,MatrixShape>& m)
      : itsLimits(m.GetLimits())
      , itsN(itsLimits.Row.size())
      , itsPtr(static_cast<SMatrix&>(m).priv_begin())
      {};

    T& operator()(index_t i,index_t j) {return itsPtr[GetSymOffset(i,j,itsLimits)];}
  private:
    const MatLimits  itsLimits;
    const index_t    itsN;
    T*               itsPtr;
  };


#undef CHECK

 private:
  friend class      Indexable<T,SMatrix,Symmetric,Real,MatrixShape>;
  friend class ArrayIndexable<T,SMatrix,Symmetric     ,MatrixShape>;
  friend class Subscriptor;

  const T* priv_begin() const;
        T* priv_begin()      ;
  void  Check () const; //Check internal consistency between limits and cow.
  static size_t  size(size_t N) {return (N*(N+1))/2;}

  size_t      itsN   ;   //NxN matrix.
#ifdef OML_USE_STDVEC
  std::vector<T> itsData;
#else
  cow_array<T> itsData;   //Copy-On-Write array for the data.
#endif
};

template <class T, class A, Data D> inline
std::ostream& operator<<(std::ostream& os,const Indexable<T,A,Symmetric,D,MatrixShape>& a)
{
  return os << SMatrix<T>(a);
}


//----------------------------------------------------------------------
//
//  Construction zone.
//
template <class T> template <class B, Data D> inline
SMatrix<T>::SMatrix(const Indexable<T,B,Symmetric,D,MatrixShape>& m)
  : MatrixBase(m.GetLimits           ())
  , itsN      (GetLimits().GetNumRows())
  , itsData   (size(itsN)           )
  {
    MatrixAssign(*this,m); //Use op[] or op(i,j) depending on D.
  }

template <class T> template <class B, Data D> inline
SMatrix<T>& SMatrix<T>::operator=(const Indexable<T,B,Symmetric,D,MatrixShape>& m)
{
  SetLimits(m.GetLimits());
  MatrixAssign(*this,m); //Use op[] or op(i,j) depending on D.
  return *this;
}

template <class T> inline size_t  SMatrix<T>::size() const
{
  return itsData.size();
}

template <class T> inline const T* SMatrix<T>::priv_begin() const
{
  return &*itsData.begin();
}

template <class T> inline T* SMatrix<T>::priv_begin()
{
  return &*itsData.begin();
}

template <class T> inline void SMatrix<T>::SetLimits(size_t N, bool preserve)
{
  SetLimits(MatLimits(N,N),preserve);
}

template <class T> inline void SMatrix<T>::SetLimits(const VecLimits& lim, bool preserve)
{
  SetLimits(MatLimits(lim,lim),preserve);
}

template <class T> inline MatLimits SMatrix<T>::GetLimits() const
{
  return MatrixBase::GetLimits();
}
//namespace std {
//template <class T> class complex;
//};
typedef std::complex<double> dcmplx;

template <class T> inline  const T& SMatrix<T>::ref(index_t i,index_t j) const
{
    index_t index=GetSymOffset(i,j,GetLimits());
    return itsData[index];
}

template <class T> inline  T SMatrix<T>::operator()(index_t i,index_t j) const
{
    index_t index=GetSymOffset(i,j,GetLimits());
    return itsData[index];
}
template <class T> inline  T& SMatrix<T>::operator()(index_t i,index_t j)
{
    assert(i<=j);
    index_t index=GetSymOffset(i,j,GetLimits());
    return itsData[index];
}

template <> inline  dcmplx SMatrix<dcmplx>::operator()(index_t i,index_t j) const
{
    index_t index=GetSymOffset(i,j,GetLimits());
    return i<=j ? itsData[index] : conj(itsData[index]);
}
template <> inline  dcmplx& SMatrix<dcmplx>::operator()(index_t i,index_t j)
{
    assert(i<=j);
    index_t index=GetSymOffset(i,j,GetLimits());
    return itsData[index];
}


inline void Fill(SMatrix<dcmplx>& a, const dcmplx& val)
{
    for (index_t i:a.rows())
    {
        a(i,i)=real(val);
        for (index_t j:a.cols(i+1))
            a(i,j)=val;
    }
}

template <class T> inline void Unit(SMatrix<T>& a)
{
    for (index_t i:a.rows())
    {
        a(i,i)=T(1);
        for (index_t j:a.cols(i+1))
            a(i,j)=T(0);
    }
}

template <class T, class A, Data D> inline
double FrobeniusNorm(const Indexable<T,A,Symmetric,D,MatrixShape>& m)
{
    double fnorm=0.0;
    for (index_t i: m.rows())
    {
        fnorm+=real(m(i,i)*conj(m(i,i)));
        for (index_t j: m.cols(i+1))
            fnorm+=real(m(i,j)*conj(m(i,j)));
    }
    return sqrt(fnorm);
}

//---------------------------------------------------------------------
//
//  SMatrix * SMatrix, returns a proxy.
//
template <class T, class A, class B> class MatrixSSOp
: public Indexable<T,MatrixSSOp<T,A,B>,Full,Abstract,MatrixShape>
{
 public:
  typedef Indexable<T,MatrixSSOp<T,A,B>,Full,Abstract,MatrixShape> IndexableT;
  typedef Ref<T,IndexableT,MatrixShape> RefT;

  MatrixSSOp(const A& a, const B& b)
    : itsA(a)
    , itsB(b)
  {
    //std::cout << "Alim="<< itsA.GetLimits() << " Blim=" << itsB.GetLimits()<< std::endl;
    assert(itsA.GetLimits().Col==itsB.GetLimits().Row);
  };
  MatrixSSOp(const MatrixSSOp& m)
    : itsA(m.itsA)
    , itsB(m.itsB)
  {}
  T operator()(index_t i, index_t j) const
  {
    T ret(0);
    for (index_t k:itsA.cols())
        ret+=itsA(i,k)*itsB(k,j);
    return ret;
  }
  size_t  size() const {return MatLimits().size();}
  MatLimits GetLimits() const {return MatLimits(itsA.GetLimits().Row,itsB.GetLimits().Col);}

 private:
  const A itsA;
  const B itsB;
};

template <class TA, class TB, class A, class B, Data DA, Data DB> inline
auto operator*(const Indexable<TA,A,Symmetric,DA,MatrixShape>& a,const Indexable<TB,B,Symmetric,DB,MatrixShape>& b)
{
  typedef typename ReturnType<TA,TB>::RetType TR;
  return MatrixSSOp<TR,typename A::RefT,typename B::RefT>(a,b);
}
template <class TA, class TB, class A, class B, Data DA, Data DB> inline
auto operator*(const Indexable<TA,A,Full,DA,MatrixShape>& a,const Indexable<TB,B,Symmetric,DB,MatrixShape>& b)
{
  typedef typename ReturnType<TA,TB>::RetType TR;
  return MatrixSSOp<TR,typename A::RefT,typename B::RefT>(a,b);
}
template <class TA, class TB, class A, class B, Data DA, Data DB> inline
auto operator*(const Indexable<TA,A,Symmetric,DA,MatrixShape>& a,const Indexable<TB,B,Full,DB,MatrixShape>& b)
{
  typedef typename ReturnType<TA,TB>::RetType TR;
  return MatrixSSOp<TR,typename A::RefT,typename B::RefT>(a,b);
}


#endif //_SMatrix_H_
