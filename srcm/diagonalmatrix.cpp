module;
#include <complex>
#include <vector>
#include <cassert>
#include <iomanip>

export module oml.DiagonalMatrix;


export import oml.VecLimits;
export import oml.MatLimits;
export import oml.MatIndex;
export import oml.FakeDouble;
import oml.Shape;
import oml.unop;
import oml.MixTypes;
import oml.vector;
import oml.Xpr;
import oml.Indexable;
import oml.ArrIndex;
import oml.CopyOnWrite;
import oml.Matrixbase;
import oml.RowCol;
import oml.Matrix; //Need some operator definitions.

//----------------------------------------------------------------------------
//
//  Diagonal matrix class.  For now this is restricted to square matrix shape.
//
export template <class T> class DiagonalMatrix
    : public MatrixBase
    , public ArrayIndexable<T,DiagonalMatrix<T>,Diagonal,MatrixShape  >
    , public Indexable<T,DiagonalMatrix<T>,Diagonal,Real,MatrixShape>
    , public TStreamableObject<DiagonalMatrix<T> >
{
public:
    typedef Indexable<T,DiagonalMatrix<T>,Diagonal,Real,MatrixShape> IndexableT;
    typedef Ref<T,IndexableT,MatrixShape> RefT;

    explicit DiagonalMatrix();
    explicit DiagonalMatrix(size_t n);
    explicit DiagonalMatrix(const VecLimits&);
    explicit DiagonalMatrix(const MatLimits&);
    explicit DiagonalMatrix(const Vector<T>&);
    DiagonalMatrix(const DiagonalMatrix&);
    DiagonalMatrix& operator=(const DiagonalMatrix& dm)
    {
        itsData=dm.itsData;
        MatrixBase::SetLimits(dm.GetLimits());
        return *this;
    }
    DiagonalMatrix& operator=(const Vector<T>& v)
    {
        itsData=v;
        MatrixBase::SetLimits(MatLimits(v.size(),v.size()));
        return *this;
    }

    template <class A>         DiagonalMatrix(const Indexable<T,A,Diagonal,Real,MatrixShape>&);
    template <class A, Data D> DiagonalMatrix(const Indexable<T,A,Diagonal,D,MatrixShape>&);

    template <class A>         DiagonalMatrix& operator=(const Indexable<T,A,Diagonal,Real,MatrixShape>&);
    template <class A, Data D> DiagonalMatrix& operator=(const Indexable<T,A,Diagonal,D,MatrixShape>&);

#ifdef OML_MOVE_OPS
  DiagonalMatrix(DiagonalMatrix&& m);
  DiagonalMatrix& operator=(DiagonalMatrix&&);
#endif

    MatLimits ReBase(const MatLimits& newLimits)
    {
        assert(newLimits.Row.Low==newLimits.Col.Low);
        itsData.ReBase(newLimits.Row);
        return MatrixBase::ReBase(newLimits);
    }
    MatLimits ReBase(int low)
    {
        itsData.ReBase(low);
        return MatrixBase::ReBase(low,low);
    }

    std::ostream& Write(std::ostream&) const;
    std::istream& Read (std::istream&)      ;
    inline friend std::ostream& operator<<(std::ostream& os,const DiagonalMatrix& a)
    {
        return os << static_cast<const TStreamableObject<DiagonalMatrix<T> >& >(a);
    }


    T  operator()(index_t,index_t) const;
    T& operator()(index_t        )      ;

//    T operator[](index_t n) const
//    {
//        return itsData[n];
//    }

    size_t    size  () const; //Required by iterable.
    MatLimits GetLimits() const;

    const Vector<T>& GetDiagonal() const {return itsData;}
          Vector<T>& GetDiagonal()       {return itsData;}

    void SetLimits(const MatLimits& lim,bool preserve=false)
    {
        assert(lim.Row==lim.Col);
        MatrixBase::SetLimits(lim);
        itsData.SetLimits(lim.Row,preserve);
    }

    void SetLimits(index_t n,bool preserve=false)
    {
        SetLimits(MatLimits(n,n),preserve);
    }
    void ReIndexRows   (const std::vector<index_t>& index);
    void ReIndexColumns(const std::vector<index_t>& index);
    void SwapRows   (index_t i,index_t j);
    void SwapColumns(index_t i,index_t j);

    DiagonalMatrix SubMatrix(const MatLimits& lim) const;
    DiagonalMatrix SubMatrix(size_t N) const
    {
        return SubMatrix(MatLimits(N,N));
    }

    DiagonalMatrix& operator*=(const T& s)
    {
        itsData*=s;
        return *this;
    }

//-----------------------------------------------------------------------------
//
//  Allows fast L-value access.  Does COW check during construction.
//
  class Subscriptor
  {
  public:
    Subscriptor(Indexable<T,DiagonalMatrix,Diagonal,Real,MatrixShape>& m)
      : itsLimits(m.GetLimits().Row)
      , itsPtr(static_cast<DiagonalMatrix&>(m).priv_begin())
      {};

    T& operator()(index_t i)
    {
        assert(itsLimits.CheckIndex(i));
        return itsPtr[itsLimits.Offset(i)];
    }

  private:
    VecLimits itsLimits;
    T*        itsPtr;
  };


private:
    friend class Indexable<T,DiagonalMatrix,Diagonal,Real,MatrixShape>;
    friend class ArrayIndexable <T,DiagonalMatrix,Diagonal     ,MatrixShape>;

//    T  operator[](index_t i) const {return priv_begin()[i];}
//    T& operator[](index_t i)       {return priv_begin()[i];}

    const T* priv_begin() const {return itsData.begin();}
          T* priv_begin()       {return itsData.begin();}

    Vector<T> itsData;   //Copy-On-Write array for the data.
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


template <class T> inline T DiagonalMatrix<T>::operator()(index_t i,index_t j) const
{
    static T zero(0);
    CHECK(i,j);
    return i==j ? itsData(i) : zero;
}

template <class T> inline T& DiagonalMatrix<T>::operator()(index_t i)
{
    CHECK(i,i);
    return itsData(i);
}

#undef CHECK

template <class T> inline
DiagonalMatrix<T>::DiagonalMatrix()
    : DiagonalMatrix<T>(VecLimits(0))
{}

template <class T> inline
DiagonalMatrix<T>::DiagonalMatrix(size_t N)
    : DiagonalMatrix<T>(VecLimits(N))
{}

template <class T> inline
DiagonalMatrix<T>::DiagonalMatrix(const VecLimits& lim)
    : MatrixBase(MatLimits(lim,lim))
    , itsData   (lim)
{}

template <class T> inline
DiagonalMatrix<T>::DiagonalMatrix(const MatLimits& lim)
    : DiagonalMatrix<T>(lim.Row)
{
    assert(lim.Row==lim.Col);
}

template <class T> inline
DiagonalMatrix<T>::DiagonalMatrix(const Vector<T>& v)
    : MatrixBase(MatLimits(v.size(),v.size()))
    , itsData   (v) //Shallow copy
{}

template <class T> inline
DiagonalMatrix<T>::DiagonalMatrix(const DiagonalMatrix& dm)
    : MatrixBase(dm.GetLimits())
    , itsData   (dm.itsData) //Shallow copy
{}

template <class T> template <class A>    inline      DiagonalMatrix<T>::
DiagonalMatrix(const Indexable<T,A,Diagonal,Real,MatrixShape>& dm)
    : DiagonalMatrix<T>(dm.GetLimits())
{
    ArrayAssign(*this,dm); //Use op[].
}

template <class T> template <class A, Data D> inline DiagonalMatrix<T>::
DiagonalMatrix(const Indexable<T,A,Diagonal,D,MatrixShape>& dm)
    : DiagonalMatrix<T>(dm.GetLimits())
{
    DiagonalAssign(*this,dm); //Use op(i,j).
}

template <class T> template <class A> inline DiagonalMatrix<T>& DiagonalMatrix<T>::
operator=(const Indexable<T,A,Diagonal,Real,MatrixShape>& dm)
{
    SetLimits(dm.GetLimits());
    ArrayAssign(*this,dm); //Use op[].
    return *this;
}
template <class T> template <class A, Data D> inline DiagonalMatrix<T>& DiagonalMatrix<T>::
operator=(const Indexable<T,A,Diagonal,D,MatrixShape>& dm)
{
    SetLimits(dm.GetLimits());
    DiagonalAssign(*this,dm); //Use op(i,j).
    return *this;
}

#ifdef OML_MOVE_OPS

template <class T> inline DiagonalMatrix<T>::DiagonalMatrix(DiagonalMatrix<T>&& m)
  : MatrixBase(m)
  , itsData   (std::move(m.itsData))
  {
//    std::cout << "DiagonalMatrix<T> move constructor m.itsData.size()=" << m.itsData.size() << std::endl;
  }

template <class T> inline DiagonalMatrix<T>& DiagonalMatrix<T>::operator=(DiagonalMatrix<T>&& m)
{
  MatrixBase::operator=(m);
  itsData=std::move(m.itsData);
//  std::cout << "DiagonalMatrix<T> move op=" << std::endl;
  return *this;
}
#endif

export template <class T, class Derived, Data D, class B, Data DB> inline
void DiagonalAssign(Indexable<T,Derived,Diagonal,D,MatrixShape>& a,const Indexable<T,B,Diagonal,DB,MatrixShape>& b)
{
#ifdef WARN_DEEP_COPY
  std::cerr << "Doing abstract MatrixAssign n=" << a.size() << std::endl;
#endif
  assert(a.GetLimits()==b.GetLimits());
  typename Derived::Subscriptor s(a);
  for (index_t i:a.rows())
      s(i)=b(i,i);
}

template <class T> inline size_t DiagonalMatrix<T>::size() const
{
    return itsData.size();
}

template <class T> inline MatLimits DiagonalMatrix<T>::GetLimits() const
{
    return MatrixBase::GetLimits();
}

export template <class T, class A, Shape S> inline T Sum(const ArrayIndexable<T,A,Diagonal,S>& a)
{
  T ret(0);
  for (const T& ai:a) ret+=ai;
  return ret;
}


//--------------------------------------------------------------------------
//
//  Matrix algebra functions  D*D, D*M M*D
//
//  All of these operators return a proxy for the operation.  So for example a transpose
//  proxy just stores a Matrix reference and op(i,j) just return ref(j,i). The multiplication
//  proxies are little more complicated, but essentially just return a Row(i)*Col(j) dot
//  product.
//----------------------------------------------------------------------
//
//  Multiplication, Matrix * Diagonal Proxy.
//
template <class T, class A, class D> class MatrixMDOp
: public Indexable<T,MatrixMDOp<T,A,D>,Full,Abstract,MatrixShape>
{
 public:
  typedef Indexable<T,MatrixMDOp<T,A,D>,Full,Abstract,MatrixShape> IndexableT;
  typedef Ref<T,IndexableT,MatrixShape> RefT;

  MatrixMDOp(const A& a, const D& d)
    : itsA(a)
    , itsD(d)
  {
    assert(itsA.GetLimits().Col==itsD.GetLimits().Row);
  };
  MatrixMDOp(const MatrixMDOp& m)
    : itsA(m.itsA)
    , itsD(m.itsD)
    {};
  T operator()(index_t i,index_t j) const
  {
    return itsA(i,j)*itsD(j,j);
  }
  MatLimits GetLimits() const {return MatLimits(itsA.GetLimits().Row,itsD.GetLimits().Col);}
  size_t    size  () const {return GetLimits().size();}

 private:
  const A itsA;
  const D itsD;
};

//----------------------------------------------------------------------
//
//  Multiplication, Diagonal * Matrix Proxy.
//
template <class T, class D, class B> class MatrixDMOp
: public Indexable<T,MatrixDMOp<T,D,B>,Full,Abstract,MatrixShape>
{
 public:
  typedef Indexable<T,MatrixDMOp<T,D,B>,Full,Abstract,MatrixShape> IndexableT;
  typedef Ref<T,IndexableT,MatrixShape> RefT;

  MatrixDMOp(const D& d, const B& b)
    : itsD(d)
    , itsB(b)
  {
    assert(itsD.GetLimits().Col==itsB.GetLimits().Row);
  };
  MatrixDMOp(const MatrixDMOp& m)
    : itsD(m.itsD)
    , itsB(m.itsB)
    {};
  T operator()(index_t i,index_t j) const
  {
    return itsD(i,i)*itsB(i,j);
  }
  MatLimits GetLimits() const {return MatLimits(itsD.GetLimits().Row,itsB.GetLimits().Col);}
  size_t    size  () const {return GetLimits().size();}

 private:
  const D itsD;
  const B itsB;
};

//----------------------------------------------------------------------
//
//  Multiplication, Diagonal * Diagonal Proxy.
//

template <class T, class D, class B> class MatrixDDOp
: public Indexable<T,MatrixDDOp<T,D,B>,Diagonal,Abstract,MatrixShape>
{
 public:
  typedef Indexable<T,MatrixDDOp<T,D,B>,Diagonal,Abstract,MatrixShape> IndexableT;
  typedef Ref<T,IndexableT,MatrixShape> RefT;

  MatrixDDOp(const D& d, const B& b)
    : itsD(d)
    , itsB(b)
  {
    assert(itsD.GetLimits().Col==itsB.GetLimits().Row);
  };
  MatrixDDOp(const MatrixDDOp& m)
    : itsD(m.itsD)
    , itsB(m.itsB)
    {};
  T operator()(index_t i,index_t j) const
  {
    return j==i ? itsD(i,i)*itsB(i,i) : T(0);
  }
  MatLimits GetLimits() const {return MatLimits(itsD.GetLimits().Row,itsB.GetLimits().Col);}
  size_t    size  () const {return GetLimits().size();}

 private:
  const D itsD;
  const B itsB;
};

//----------------------------------------------------------------------
//
//  Multiplication, Diagonal * Vector Proxy.
//

template <class T, class D, class V> class MatrixDVOp
: public Indexable<T,MatrixDVOp<T,D,V>,Diagonal,Abstract,VectorShape>
{
 public:
  typedef Indexable<T,MatrixDVOp<T,D,V>,Diagonal,Abstract,VectorShape> IndexableT;
  typedef Ref<T,IndexableT,VectorShape> RefT;

  MatrixDVOp(const D& d, const V& v)
    : itsD(d)
    , itsV(v)
  {
    assert(itsD.GetLimits().Col==itsV.GetLimits());
  };
  MatrixDVOp(const MatrixDVOp& m)
    : itsD(m.itsD)
    , itsV(m.itsV)
    {};
  T operator()(index_t i) const
  {
    return itsD(i,i)*itsV(i);
  }
  VecLimits GetLimits() const {return itsD.GetLimits().Row;}
  size_t    size  () const {return itsD.size();}

 private:
  const D itsD;
  const V itsV;
};

//----------------------------------------------------------------------
//
//  Multiplication, Vector * Diagonal Proxy.
//
template <class T, class V, class D> class MatrixVDOp
: public Indexable<T,MatrixVDOp<T,V,D>,Full,Abstract,VectorShape>
{
 public:
  typedef Indexable<T,MatrixVDOp<T,V,D>,Full,Abstract,VectorShape> IndexableT;
  typedef Ref<T,IndexableT,VectorShape> RefT;

  MatrixVDOp(const V& v, const D& d)
    : itsV(v)
    , itsD(d)
  {
    assert(itsV.GetLimits()==itsD.GetLimits().Row);
  };
  MatrixVDOp(const MatrixVDOp& m)
    : itsV(m.itsV)
    , itsD(m.itsD)
    {};
  T operator()(index_t i) const
  {
    return itsV(i)*itsD(i,i);
  }
  VecLimits GetLimits() const {return itsD.GetLimits().Col;}
  size_t    size  () const {return itsD.size();}

 private:
  const V itsV;
  const D itsD;
};

//---------------------------------------------------------------------
//
//  Matrix * DiagonalMatrix, returns a proxy.
//
export template <class TA, class TB, class A, class B, Data DA, Data DB> inline
auto operator*(const Indexable<TA,A,Full,DA,MatrixShape>& a,const Indexable<TB,B,Diagonal,DB,MatrixShape>& b)
{
  typedef typename ReturnType<TA,TB>::RetType TR;
  return MatrixMDOp<TR,typename A::RefT,typename B::RefT>(a,b);
}

//------------------------------------------------------------
//
//  DiagonalMatrix * Matrix, returns a proxy.
//
export template <class TA,class TB, class A, class B, Data DA, Data DB> inline
auto operator*(const Indexable<TA,A,Diagonal,DA,MatrixShape>& a,const Indexable<TB,B,Full,DB,MatrixShape>& b)
{
  typedef typename ReturnType<TA,TB>::RetType TR;
  return MatrixDMOp<TR,typename A::RefT,typename B::RefT>(a,b);
}

//------------------------------------------------------------
//
//  DiagonalMatrix * DiagonalMatrix, returns a proxy.
//
export template <class TA,class TB, class A, class B, Data DA, Data DB> inline
auto operator*(const Indexable<TA,A,Diagonal,DA,MatrixShape>& a,const Indexable<TB,B,Diagonal,DB,MatrixShape>& b)
{
  typedef typename ReturnType<TA,TB>::RetType TR;
  return MatrixDDOp<TR,typename A::RefT,typename B::RefT>(a,b);
}

export template <class TA,class TB, class A, class B, Data DA, Data DB> inline
auto operator*(const Indexable<TA,A,Diagonal,DA,MatrixShape>& a,const Indexable<TB,B,Full,DB,VectorShape>& b)
{
  typedef typename ReturnType<TA,TB>::RetType TR;
  return MatrixDVOp<TR,typename A::RefT,typename B::RefT>(a,b);
}

export template <class TA,class TB, class A, class B, Data DA, Data DB> inline
auto operator*(const Indexable<TA,A,Full,DA,VectorShape>& a,const Indexable<TB,B,Diagonal,DB,MatrixShape>& b)
{
  typedef typename ReturnType<TA,TB>::RetType TR;
  return MatrixVDOp<TR,typename A::RefT,typename B::RefT>(a,b);
}

//-----------------------------------------------------------------------------
//
//  IO
//
template <class T> std::ostream& DiagonalMatrix<T>::Write(std::ostream& os) const
{
  assert(os);
  if (!this->Pretty())
  {
    os << GetLimits();
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
        os << std::setw(wid) << std::setprecision(prec) << (*this)(i,j) << " ";
      os << "]" << std::endl;
    }
  }
  return os;
}

template <class T> std::istream& DiagonalMatrix<T>::Read(std::istream& is)
{
  assert(is);
  MatLimits lim;
  is >> lim;
  if (size()==0)
    SetLimits(lim.size());
   else
    assert(size()==0);

  ::Read(is,*this);
  assert(is);
  return is;
}

template <class T> void DiagonalMatrix<T>::ReIndexRows(const std::vector<index_t>& index)
{
    GetDiagonal().ReIndex(index);
}

template <class T> void DiagonalMatrix<T>::ReIndexColumns(const std::vector<index_t>& index)
{
    GetDiagonal().ReIndex(index);
}

template <class T> void DiagonalMatrix<T>::SwapRows(index_t i,index_t j)
{
     T temp=(*this)(i,i);
     (*this)(i)=(*this)(j,j);
     (*this)(j)=temp;
}

template <class T> void DiagonalMatrix<T>::SwapColumns(index_t i,index_t j)
{
     T temp=(*this)(i,i);
     (*this)(i)=(*this)(j,j);
     (*this)(j)=temp;
}


template <class T> DiagonalMatrix<T> DiagonalMatrix<T>::SubMatrix(const MatLimits& lim) const
{
	assert(lim.Row.Low >=GetLimits().Row.Low );
	assert(lim.Row.High<=GetLimits().Row.High);
	DiagonalMatrix<T> dest(lim);
	for (index_t i:dest.rows())
        dest(i)=(*this)(i,i);
    return dest;
}


