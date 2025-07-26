module;
#include <complex>
#include <vector>
#include <cassert>
#include <iomanip>

export module oml.SMatrix;


export import oml.VecLimits;
export import oml.MatLimits;
export import oml.MatIndex;
export import oml.FakeDouble;
import oml.Shape;
import oml.unop;
import oml.MixTypes;
import oml.Vector;
import oml.Xpr;
import oml.Indexable;
import oml.ArrIndex;
import oml.CopyOnWrite;
import oml.Matrixbase;
import oml.RowCol;
import oml.Matrix; //Need some operator definitions.
export import oml.SMatrixIndex; //Just to get Sum?

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
export template <class T> class SMatrix
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
  inline friend std::ostream& operator<<(std::ostream& os,const SMatrix& a)
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

export template <class T, class A, Data D> inline
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


export inline void Fill(SMatrix<dcmplx>& a, const dcmplx& val)
{
    for (index_t i:a.rows())
    {
        a(i,i)=real(val);
        for (index_t j:a.cols(i+1))
            a(i,j)=val;
    }
}

export template <class T> inline void Unit(SMatrix<T>& a)
{
    for (index_t i:a.rows())
    {
        a(i,i)=T(1);
        for (index_t j:a.cols(i+1))
            a(i,j)=T(0);
    }
}

export template <class T, class A, Data D> inline
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

export template <class TA, class TB, class A, class B, Data DA, Data DB> inline
auto operator*(const Indexable<TA,A,Symmetric,DA,MatrixShape>& a,const Indexable<TB,B,Symmetric,DB,MatrixShape>& b)
{
  typedef typename ReturnType<TA,TB>::RetType TR;
  return MatrixSSOp<TR,typename A::RefT,typename B::RefT>(a,b);
}
export template <class TA, class TB, class A, class B, Data DA, Data DB> inline
auto operator*(const Indexable<TA,A,Full,DA,MatrixShape>& a,const Indexable<TB,B,Symmetric,DB,MatrixShape>& b)
{
  typedef typename ReturnType<TA,TB>::RetType TR;
  return MatrixSSOp<TR,typename A::RefT,typename B::RefT>(a,b);
}
export template <class TA, class TB, class A, class B, Data DA, Data DB> inline
auto operator*(const Indexable<TA,A,Symmetric,DA,MatrixShape>& a,const Indexable<TB,B,Full,DB,MatrixShape>& b)
{
  typedef typename ReturnType<TA,TB>::RetType TR;
  return MatrixSSOp<TR,typename A::RefT,typename B::RefT>(a,b);
}


#ifdef DEBUG
#define CHECK Check()
#else
#define CHECK
#endif

//-----------------------------------------------------------------------------
//
//  Constructors.
//
template <class T> SMatrix<T>::SMatrix()
  : itsN   (0)
  , itsData(0)
  {}

template <class T> SMatrix<T>::SMatrix(index_t N)
  : MatrixBase(MatLimits(N,N))
  , itsN      (N)
  , itsData   (size(itsN))
  {
    CHECK;
  }

template <class T> SMatrix<T>::SMatrix(const VecLimits& r,const VecLimits& c)
  : MatrixBase(MatLimits(r,c))
  , itsN      (GetLimits().GetNumRows() )
  , itsData   (size(itsN))
  {
    CHECK;
  }

template <class T> SMatrix<T>::SMatrix(const MatLimits& lim)
  : MatrixBase(lim              )
  , itsN      (GetLimits().GetNumRows() )
  , itsData   (size(itsN))
  {
    CHECK;
  }

template <class T> SMatrix<T>::SMatrix(const SMatrix<T>& m)
  : MatrixBase(m)
  , itsN      (m.itsN)
  , itsData   (m.itsData)
  {
    CHECK;
  }
	
template <class T> SMatrix<T>::SMatrix(const Matrix<T>& m, double tol)
  : MatrixBase(m.GetLimits())
  , itsN      (GetLimits().GetNumRows() )
  , itsData   (size(itsN))
  {
    CHECK;
    for (index_t i:this->rows())
        for (index_t j:this->cols(i))
        {
            assert(Max(fabs(m(i,j)-conj(m(j,i))))<=tol);
            (*this)(i,j)=m(i,j);
        }
  }


template <class T> SMatrix<T>& SMatrix<T>::operator=(const SMatrix<T>& m)
{
  MatrixBase::operator=(m);
  itsData=m.itsData;
  itsN=m.itsN;
  CHECK;
  return *this;
}

#undef CHECK


//-----------------------------------------------------------------------------
//
//  Internal consistency check.
//
template <class T> void SMatrix<T>::Check() const
{
  assert(GetLimits().Check());
  assert(GetRowLimits()==GetColLimits());
  assert(itsN==GetNumRows());
  assert(size(itsN)==itsData.size());
}


//-----------------------------------------------------------------------------
//
//  IO
//
template <class T> std::ostream& SMatrix<T>::Write(std::ostream& os) const
{
  assert(os);
  if (!this->Pretty())
  {
    os << GetLimits();
    if (this->Binary())
      BinaryWrite(itsN,os);
    else
      os << itsN << " ";
    ::Write(os,*this);
  }
  else
  {
    int prec=os.precision();
    int wid =os.width();
    os << std::setw(0) << GetLimits() << std::endl;
    for (index_t i=GetRowLow();i<=GetRowHigh();i++)
    {
      os << "[ ";
      for (index_t j=GetColLow();j<=GetColHigh();j++)
      {
//        if (i<=j)
	  os << std::setw(wid) << std::setprecision(prec) << (*this)(i,j) << " ";
//	else
//	  for (int n=0;n<=wid;n++) os << " ";
      }
      os << "]" << std::endl;
    }
  }
  return os;
}

template <class T> std::istream& SMatrix<T>::Read(std::istream& is)
{
  assert(is);
  MatLimits lim;
  is >> lim;

  if (this->Binary())
    BinaryRead(itsN,is);
  else
    is >> itsN;

  if (size()==0)
    SetLimits(lim);
   else
    assert(GetLimits()==lim);

  ::Read(is,*this);
#if DEBUG
	Check();
#endif
  assert(is);
  return is;
}

//-----------------------------------------------------------------------------
//
//  Changing size and/or limits.
//
template <class T> void SMatrix<T>::SetLimits(const MatLimits& theLimits, bool preserve)
{
#ifdef DEBUG
  theLimits.Check();
#endif
  if (GetLimits()!=theLimits)
  {
    if (preserve)
    {
      SMatrix<T> dest(theLimits);

      index_t newrowlow=theLimits.Row.Low, newrowhigh=theLimits.Row.High;
      //index_t newcollow=theLimits.Col.Low;
      index_t newcolhigh=theLimits.Col.High;
      index_t rowlow =Max(newrowlow ,GetLimits().Row.Low );   //Limits of old and new
      index_t rowhigh=Min(newrowhigh,GetLimits().Row.High);   //data overlap.
//      index_t collow =Max(newcollow ,GetLimits().Col.Low );   //Limits of old and new
      index_t colhigh=Min(newcolhigh,GetLimits().Col.High);   //data overlap.

      Subscriptor sdest (dest);         //Destination subscriptor.

      for (index_t i=rowlow;i<=rowhigh;i++)
        for (index_t j=i;j<=colhigh;j++)
          sdest(i,j)=(*this)(i,j); //Transfer any overlaping data.

			(*this)=dest;
		}
     else
    {
			(*this)=SMatrix<T>(theLimits);
    }
  }
}


template <class T> void SMatrix<T>::ReIndexRows(const std::vector<index_t>& index)
{
  assert(GetLimits().GetNumRows()==index.size());

  typename std::vector<index_t>::const_iterator i=index.begin();
  SMatrix<T> dest(GetLimits());
  Subscriptor sdest(dest);

  for (int row=GetLimits().Row.Low;row<=GetLimits().Row.High;row++,i++)
  for (int col=GetLimits().Col.Low;col<=GetLimits().Col.High;col++)
  sdest(row,col)=(*this)(*i+GetLimits().Row.Low,col);

  *this=dest;
}

template <class T> void SMatrix<T>::ReIndexColumns(const std::vector<index_t>& index)
{
  assert(GetLimits().GetNumCols()==index.size());

  typename std::vector<index_t>::const_iterator i=index.begin();
  SMatrix<T> dest(GetLimits());
  Subscriptor sdest(dest);

  for (int col=GetLimits().Col.Low;col<=GetLimits().Col.High;col++,i++)
    for (int row=GetLimits().Row.Low;row<=GetLimits().Row.High;row++)
      sdest(row,col)=(*this)(row,*i+GetLimits().Col.Low);

  *this=dest;
}
//-----------------------------------------------------------------------------
//
//  Swapping
//

template <class T> void SMatrix<T>::SwapRows(index_t i,index_t j)
{
  Subscriptor s(*this);
  for (index_t c : this->cols())
  {
     T temp=s(i,c);
     s(i,c)=s(j,c);
     s(j,c)=temp;
  }
}

template <class T> void SMatrix<T>::SwapColumns(index_t i,index_t j)
{
  Subscriptor s(*this);
  for (index_t r : this->rows())
  {
     T temp=s(r,i);
     s(r,i)=s(r,j);
     s(r,j)=temp;
  }
}
template <class T> SMatrix<T> SMatrix<T>::SubMatrix(const MatLimits& lim) const
{
    SMatrix<T> dest(lim);
	assert(dest.GetLimits().Row.Low >=GetLimits().Row.Low );
	assert(dest.GetLimits().Col.Low >=GetLimits().Col.Low );
	assert(dest.GetLimits().Row.High<=GetLimits().Row.High);
	assert(dest.GetLimits().Col.High<=GetLimits().Col.High);
	Subscriptor s(dest);
	index_t rh=dest.GetLimits().Row.High;
	index_t ch=dest.GetLimits().Col.High;
	for (index_t i=dest.GetLimits().Row.Low;i<=rh;i++)
		for (index_t j=i;j<=ch;j++)
			s(i,j)=(*this)(i,j);
    return dest;
}

