// File: matrixm.cpp  Make a module for Matrix<T>
module;
#include <cstddef>
#include <iostream>
#include <cassert>
#include <vector>
#include <cmath>
#include <iomanip>

export module oml.matrix;
import oml.vector;

//---------------------------------------------------------------------------
/*! \class MatLimits matlimit.h oml/matlimit.h
  \brief Encapsulate all effects of matrix subscript bounds in this class.

  Structure for matrix dimension limits.  All limits are checked in the
  VecLimits constructors.

  \nosubgrouping
*/
export class MatLimits
{
 public:
  /*! \name Constructors*/
  //@{
  MatLimits(                                 );//!<Limits for a null Matrix.
  MatLimits(index_t,index_t,index_t,index_t  );//!<Construct from lower and upper bounds.
  MatLimits(size_t,size_t                    );//!<Construct from row and column size, use default lower bound.
  MatLimits(const VecLimits&,const VecLimits&);//!<Construct from lower and upper bounds.
  //@}
 ~MatLimits();

  void ReBase(int rowLow,int colLow);

  static size_t  size(size_t,size_t                  );
  static size_t  size(index_t,index_t,index_t,index_t  );
  static size_t  size(const VecLimits&,const VecLimits&);

  //! Returns number of rows.
  size_t  GetNumRows        (       ) const;
  //! Returns number of columns.
  size_t  GetNumCols        (       ) const;
  //! Returns number of elements.
  size_t  size              (       ) const;

  index_t Offset    (index_t,index_t) const;
  bool    Check     (               ) const;
  bool    CheckIndex(index_t,index_t) const;

  //! Comparison.
  bool operator==(const MatLimits&) const;
  //! Comparison.
  bool operator!=(const MatLimits&) const;

  /*! \name IO
    See StreamableObject for details on output formats. You must include matrix_io.h to get
    definitions.  op>> and op<< are also defined.
  */
  //@{
  std::ostream& Write(std::ostream&) const;//!<Write to stream.
  std::istream& Read (std::istream&)      ;//!<Read from stream.
  //@}

  VecLimits Row; //!<Row    index limits.
  VecLimits Col; //!<Column index limits.
};


//------------------------------------------------------------------------
//
//  Constructors.
//
inline MatLimits::MatLimits()
  : Row()
  , Col()
  {}

inline MatLimits::MatLimits(index_t rowLow, index_t rowHigh,index_t colLow, index_t colHigh)
  : Row(rowLow,rowHigh)
  , Col(colLow,colHigh)
  {}

inline MatLimits::MatLimits(size_t rowSize, size_t colSize)
  : Row(rowSize)
  , Col(colSize)
  {}

inline MatLimits::MatLimits(const VecLimits& rowLimits, const VecLimits& colLimits)
  : Row(rowLimits)
  , Col(colLimits)
  {}

inline MatLimits::~MatLimits() {}

inline void MatLimits::ReBase(int rowLow,int colLow)
{
    Row.ReBase(rowLow);
    Col.ReBase(colLow);
}

//------------------------------------------------------------------------
//
//  Get size static functions.
//
inline size_t  MatLimits::size(size_t numRows, size_t numCols)
{
  return numRows*numCols;
}
inline size_t  MatLimits::size(const VecLimits& row, const VecLimits& col)
{
  return row.size()*col.size();
}

inline size_t  MatLimits::size(index_t rowLow, index_t rowHigh,index_t colLow, index_t colHigh)
{
  return size(VecLimits(rowLow,rowHigh),VecLimits(colLow,colHigh) );
}

inline size_t  MatLimits::GetNumRows() const
{
  return Row.size();
}

inline size_t  MatLimits::GetNumCols() const
{
  return Col.size();
}

inline size_t  MatLimits::size() const
{
  return size(GetNumRows(),GetNumCols());
}

//------------------------------------------------------------------------
//
//  Overloaded ==
//
inline bool MatLimits::operator==(const MatLimits& lim) const
{
  return (Row==lim.Row)&&(Col==lim.Col);
}

inline bool MatLimits::operator!=(const MatLimits& l) const
{
  return !((*this)==l);
}

//-------------------------------------------------------------------
//
//  Other inline functions for the MatLimits class.
//
inline bool MatLimits::Check() const
{
  return Row.Check()&&Col.Check();
}

inline index_t MatLimits::Offset(index_t i,index_t j) const
{
  return Row.Offset(i)+GetNumRows()*Col.Offset(j);
}

inline std::ostream& operator<<(std::ostream& os,const MatLimits& lim)
{
  return lim.Write(os);
}

inline std::istream& operator>>(std::istream& is, MatLimits& lim)
{
  return lim.Read (is);
}

// Use this to get limits of a tensor product
inline MatLimits operator*(const MatLimits& a, const MatLimits& b)
{
    return MatLimits(a.Row*b.Row,a.Col*b.Col);
}

// matrixbase.h

//----------------------------------------------------------------------------
/*! \class MatrixBase matrixbase.h oml/matrixbase.h
  \brief Provids all sorts of information on the row and column limits.
*/
class MatrixBase
{
 public:
  MatrixBase(                    );
  MatrixBase(const MatLimits& lim);
  // use dault op= and copy contructors
  /*! \name Limits info
  */
  //@{
  //! Extract info on Matrix size/shape.
  MatLimits GetLimits   () const;
  VecLimits GetRowLimits() const;
  VecLimits GetColLimits() const;
  size_t    GetNumRows  () const;
  size_t    GetNumCols  () const;
  size_t    size        () const;
  index_t   GetRowLow   () const;
  index_t   GetRowHigh  () const;
  index_t   GetColLow   () const;
  index_t   GetColHigh  () const;
  bool      IsSquare    () const;
  //@}
  void SetLimits(const MatLimits&, bool preserve=false);
  MatLimits ReBase(int rowLow,int colLow);
  MatLimits ReBase(const MatLimits&);

 private:
  MatLimits    itsLimits; //Manages the upper and lower std::vector limits.
};

inline MatrixBase::MatrixBase(                    ) : itsLimits(   ) {}
inline MatrixBase::MatrixBase(const MatLimits& lim) : itsLimits(lim) {}

inline void MatrixBase::SetLimits(const MatLimits& lim,bool /*preserve*/)
{
  itsLimits=lim;
}
inline MatLimits MatrixBase::ReBase(int rowLow,int colLow)
{
    MatLimits oldLimits=itsLimits;
    itsLimits.ReBase(rowLow,colLow);
    return oldLimits;
}
inline MatLimits MatrixBase::ReBase(const MatLimits& newLimits)
{
    MatLimits oldLimits=itsLimits;
    itsLimits.ReBase(newLimits.Row.Low,newLimits.Col.Low);
    return oldLimits;
}


inline  MatLimits MatrixBase::GetLimits() const
{
  return itsLimits;
}

inline size_t  MatrixBase::GetNumRows() const
{
  return itsLimits.Row.size();
}

inline size_t  MatrixBase::GetNumCols() const
{
  return itsLimits.Col.size();
}

inline size_t  MatrixBase::size() const
{
  return itsLimits.size();
}

inline bool MatrixBase::IsSquare() const
{
    return GetNumRows()==GetNumCols();
}


inline VecLimits MatrixBase::GetRowLimits() const
{
  return itsLimits.Row;
}

inline VecLimits MatrixBase::GetColLimits() const
{
  return itsLimits.Col;
}

inline index_t MatrixBase::GetRowLow() const
{
  return itsLimits.Row.Low;
}

inline index_t MatrixBase::GetRowHigh() const
{
  return itsLimits.Row.High;
}

inline index_t MatrixBase::GetColLow() const
{
  return itsLimits.Col.Low;
}

inline index_t MatrixBase::GetColHigh() const
{
  return itsLimits.Col.High;
}

// binop.h

template <class TR, class TA, class TB, class A, class B, Shape S> class XprBinary
{};

template <class TR, class TA, class TB, class A, class B> class XprBinary<TR,TA,TB,A,B,MatrixShape>
{
 public:
   XprBinary(const A& a,const B& b,TR(*f)(const TA&,const TB&)) : itsA(a), itsB(b), itsF(f) {};
  ~XprBinary() {};

  TR         operator[](index_t n          ) const {return itsF(itsA[n],itsB[n]);}
  TR         operator()(index_t i,index_t j) const {return itsF(itsA(i,j),itsB(i,j));}
  size_t     size      (                   ) const {return itsA.size();}
  MatLimits  GetLimits (                   ) const {return itsA.GetLimits();}
 private:
   A itsA;
   B itsB;
   TR(*itsF)(const TA&,const TB&);
};

// xpr.h
template <class T, class TR, class A, Shape S> class XprUnary
{};

template <class T, class TR, class A> class XprUnary<T,TR,A,MatrixShape>
{
 public:
  XprUnary(const A& a,TR(*f)(const T&)) : itsA(a), itsF(f) {};
  ~XprUnary() {};

  TR        operator()(index_t i,index_t j) const {return itsF(itsA(i,j));}
  size_t    size      (                   ) const {return itsA.size();}
  MatLimits GetLimits (                   ) const {return itsA.GetLimits();}
 private:
   A itsA;
   TR(*itsF)(const T&);
};


template <class T, class R> class Ref<T,R,MatrixShape>
: public IndexableBase<Ref<T,R,MatrixShape>,MatrixShape> //Get index iterators
{
 public:
  Ref(const R& r) : itsRef(r) {};
  T         operator()(index_t i,index_t j) const {return itsRef(i,j);}
  size_t    size      (                   ) const {return itsRef.size();}
  MatLimits GetLimits (                   ) const {return itsRef.GetLimits();}
 private:
  const R& itsRef;
};

template <class T, class A, Shape> class Val {};

template <class T, class A> class Val<T,A,MatrixShape>
{
 public:
  explicit Val(const T& v, const A& a) : itsVal(v), itsA(a),itsLimits(itsA.GetLimits()) {};
  T         operator()(index_t,index_t) const {return itsVal;}
  index_t   size      (               ) const {return itsLimits.size();} //itsA may be gone
  MatLimits GetLimits (               ) const {return itsLimits;} //Apparently itsA.GetLimits() fails at runtime
 private:
  T         itsVal;
  const A&  itsA;
  MatLimits itsLimits;
};

template <class T, class Expression,Store M,Data D,Shape S> class Xpr;

template <class T, class Expression,Store M,Data D> class Xpr<T,Expression,M,D,MatrixShape>
: public Indexable<T,Xpr<T,Expression,M,D,MatrixShape>,M,Abstract,MatrixShape>
{
 public:
  typedef Indexable<T,Xpr<T,Expression,M,D,MatrixShape>,M,Abstract,MatrixShape> IndexableT;
  typedef Ref<T,IndexableT,MatrixShape> RefT;

  Xpr(Expression e) : itsExp(e) {};
  Xpr(const Xpr& x) : itsExp(x.itsExp) {};
  ~Xpr() {};

  T         operator()(index_t i,index_t j) const {return itsExp(i,j);}
  size_t    size      (                   ) const {return itsExp.size();}
  MatLimits GetLimits (                   ) const {return itsExp.GetLimits();}

 private:
  Expression itsExp;
};
// matindex.h

//
// Template specialization provides index iterators for vector shape
//
template <class Derived> class IndexableBase<Derived,MatrixShape>
{
    public:
    //
//  Support range based iteration for rows and columns so client code and do
//     for (index_t i:A.rows())
//        for (index_t j:A.cols())
//          {do something with A(i,j)
//
//  For something like a symmtric matrix do
//     for (index_t i:A.rows())
//        for (index_t j:A.cols(i)) //start at i
//          {do something with A(i,j)
//
//
  class index_iterator
  {
    public:
        index_iterator(index_t i) : current{i} {};
        index_iterator operator++()
        {
            index_iterator ret(*this);
            current++;
            return ret;
        }
        const index_t operator*() const {return current;}
              index_t operator*()       {return current;}
        friend bool operator!=(const index_iterator& a, const index_iterator& b)
        {
            return a.current!=b.current;
        }
    private:
        index_t current;
  };

    class iterator_proxy
    {
    public:
        iterator_proxy(const VecLimits& lim) : low(lim.Low) , high(lim.High)
        {
            if (high<low) high=low-1; //Set to terminate with no iterations
        };
        iterator_proxy(index_t l, index_t h) : low(l), high(h)
        {
            if (high<low) high=low-1; //Set to terminate with no iterations
        };
        index_iterator begin() const {return low;}
        index_iterator end  () const {return high+1;}
    private:
        index_t low;
        index_t high;
    };

  iterator_proxy rows(         ) const {return iterator_proxy(static_cast<const Derived*>(this)->GetLimits().Row);}
  iterator_proxy cols(         ) const {return iterator_proxy(static_cast<const Derived*>(this)->GetLimits().Col);}
  iterator_proxy rows(index_t j) const {return iterator_proxy(j,static_cast<const Derived*>(this)->GetLimits().Row.High);}
  iterator_proxy cols(index_t i) const {return iterator_proxy(i,static_cast<const Derived*>(this)->GetLimits().Col.High);}
  iterator_proxy array_indices (         ) const {return iterator_proxy(0,static_cast<const Derived*>(this)->size()-1);}
};

template <class T, class A, class B> class MatrixMMOp;
//-------------------------------------------------
//
//  template specialization for Matrices's.
//
template <class T, class Derived, Data D> class Indexable<T,Derived,Full,D,MatrixShape>
 : public IndexableBase<Derived,MatrixShape>
{
 public:
  explicit Indexable() {};
  ~Indexable() {};

  T  operator()(index_t i,index_t j) const {return static_cast<const Derived*>(this)->operator()(i,j);}

  size_t    size  () const {return static_cast<const Derived*>(this)->size();}
  MatLimits GetLimits() const {return static_cast<const Derived*>(this)->GetLimits();}

 private:
  Indexable& operator=(const Indexable&);
  Indexable(const Indexable&);
};
//-------------------------------------------------
//
//  template specialization for abstract Matrices's.
//
template <class T, class Derived> class Indexable<T,Derived,Full,Abstract,MatrixShape>
 : public IndexableBase<Derived,MatrixShape>
{
 public:
  explicit Indexable() {};
  ~Indexable() {};

  T operator()(index_t i,index_t j) const {return static_cast<const Derived*>(this)->operator()(i,j);}

  size_t    size  () const {return static_cast<const Derived*>(this)->size();}
  MatLimits GetLimits() const {return static_cast<const Derived*>(this)->GetLimits();}

 private:
  Indexable& operator=(const Indexable&);
  Indexable(const Indexable&);
};

//-------------------------------------------------
//
//  template specialization for diagonal Matrices's.
//
template <class T, class Derived> class Indexable<T,Derived,Diagonal,Real,MatrixShape>
 : public IndexableBase<Derived,MatrixShape>
{
 public:
  explicit Indexable() {};
  ~Indexable() {};

  T  operator()(index_t i,index_t j) const {return static_cast<const Derived*>(this)->operator()(i,j);}
  T& operator()(index_t i          )       {return static_cast<      Derived*>(this)->operator()(i);}

  size_t    size  () const {return static_cast<const Derived*>(this)->size();}
  MatLimits GetLimits() const {return static_cast<const Derived*>(this)->GetLimits();}

 private:
  Indexable& operator=(const Indexable&);
  Indexable(const Indexable&);
};

//-------------------------------------------------
//
//  template specialization for abstract Diagonal Matrices's.
//
template <class T, class Derived> class Indexable<T,Derived,Diagonal,Abstract,MatrixShape>
 : public IndexableBase<Derived,MatrixShape>
{
 public:
  explicit Indexable() {};
  ~Indexable() {};

  T operator()(index_t i,index_t j) const {return static_cast<const Derived*>(this)->operator()(i,j);}

  size_t    size  () const {return static_cast<const Derived*>(this)->size();}
  MatLimits GetLimits() const {return static_cast<const Derived*>(this)->GetLimits();}

 private:
  Indexable& operator=(const Indexable&);
  Indexable(const Indexable&);
};

//
//  Create assign functions
//
template <class T, class Derived, Data D, class B, Store MB, Data DB> inline
void MatrixAssign(Indexable<T,Derived,Full,D,MatrixShape>& a,const Indexable<T,B,MB,DB,MatrixShape>& b)
{
#ifdef WARN_DEEP_COPY
  std::cerr << "Doing abstract MatrixAssign n=" << a.size() << std::endl;
#endif
  assert(a.GetLimits()==b.GetLimits());
  typename Derived::Subscriptor s(a);
  #pragma omp parallel for collapse(2)
  for (index_t i=a.GetLimits().Row.Low;i<=a.GetLimits().Row.High;i++)
    for (index_t j=a.GetLimits().Col.Low;j<=a.GetLimits().Col.High;j++)
      s(i,j)=b(i,j);
}

template <class T, class Derived, Data D, class B, Data DB> inline
void MatrixAssign(Indexable<T,Derived,Symmetric,D,MatrixShape>& a,const Indexable<T,B,Symmetric,DB,MatrixShape>& b)
{
#ifdef WARN_DEEP_COPY
  std::cerr << "Doing abstract MatrixAssign n=" << a.size() << std::endl;
#endif
  assert(a.GetLimits()==b.GetLimits());
  typename Derived::Subscriptor s(a);
  #pragma omp parallel for collapse(2)
  for (index_t i=a.GetLimits().Row.Low;i<=a.GetLimits().Row.High;i++)
    for (index_t j=i;j<=a.GetLimits().Col.High;j++)
      s(i,j)=b(i,j);
}

#define OP(NAME,OP) \
template <class T, class Derived,Store M,Data D> inline \
Derived& Matrix##NAME (Indexable<T,Derived,M,D,MatrixShape>& a,const T& scalar)\
{\
  typename Derived::Subscriptor s(a); \
  for (index_t i:a.rows())\
    for (index_t j:a.cols())\
      s(i,j) OP##=scalar;\
  return static_cast<Derived&>(a);\
}\
template <class T, class Derived,Data D> inline \
Derived& Matrix##NAME (Indexable<T,Derived,Symmetric,D,MatrixShape>& a,const T& scalar)\
{\
 typename Derived::Subscriptor s(a); \
 for (index_t i:a.rows())\
    for (index_t j:a.cols())\
      s(i,j) OP##=scalar;\
  return static_cast<Derived&>(a);\
}\
template <class T,class Derived,Store M,Data D,class B,Data DB> inline \
Derived& Matrix##NAME (Indexable<T,Derived,M,D,MatrixShape>& a,\
                  const Indexable<T,B,M,DB,MatrixShape>& b)\
{\
    typename Derived::Subscriptor s(a); \
    for (index_t i:a.rows())\
        for (index_t j:a.cols())\
            s(i,j) OP##=b(i,j);\
  return static_cast<Derived&>(a);\
}\
template <class T,class Derived,Data D,class B,Data DB> inline \
Derived& Matrix##NAME (Indexable<T,Derived,Symmetric,D,MatrixShape>& a,\
                  const Indexable<T,B,Symmetric,DB,MatrixShape>& b)\
{\
  assert(a.GetLimits()==b.GetLimits());\
  typename Derived::Subscriptor s(a); \
  for (index_t i:a.rows())\
    for (index_t j:a.cols(i))\
      s(i,j) OP##=b(i,j);\
  return static_cast<Derived&>(a);\
}

OP(Add,+)
OP(Sub,-)
OP(Mul,*)
OP(Div,/)

#undef OP


export template <class T> class Matrix;

template <class T, class A, Store M, Data D> inline
std::ostream& operator<<(std::ostream& os,const Indexable<T,A,M,D,MatrixShape>& a)
{
  return os << Matrix<T>(a);
}

template <class T, class A> inline T Sum(const Indexable<T,A,Full,Abstract,MatrixShape>& a)
{
	T ret(0);
    for (index_t i:a.rows())
        for (index_t j:a.cols())
      ret+=a(i,j);
	return ret;
}

template <class T, class A> inline T Sum(const Indexable<T,A,Diagonal,Abstract,MatrixShape>& a)
{
  T ret(0);
  for (index_t i: a.rows()) ret+=a(i,i);
  return ret;
}


template <class T, class A, class Op, Store M, Data D,Shape S> class MinMax;

template <class T, class A, class Op, Store M> class MinMax<T,A,Op,M,Abstract,MatrixShape>
{
public:
    static T apply(const Indexable<T,A,M,Abstract,MatrixShape>& a)
    {
        int rl=a.GetLimits().Row.Low;
        int cl=a.GetLimits().Col.Low;
        T ret=a.size()>0 ? a(rl,cl) : T(0); // Don't try and read a[0] if there is no data in a!
        for (index_t i:a.rows())
        for (index_t j:a.cols())
        {
            T ai=a(i,j);
            if (Op::apply(ai,ret)) ret=ai;
        }
        return ret;
    }
};

template <class T, class A, Store M> inline T Min(const Indexable<T,A,M,Abstract,MatrixShape>& a)
{
	return MinMax<T,A,OpLT<T>,M,Abstract,MatrixShape>::apply(a);
}

template <class T, class A, Store M> inline T Max(const Indexable<T,A,M,Abstract,MatrixShape>& a)
{
	return MinMax<T,A,OpGT<T>,M,Abstract,MatrixShape>::apply(a);
}


// rowcol.h
//---------------------------------------------------------------------
//
//  Matrix row proxy class.
//
template <class T, class A, Store M, Data D> class MatrixRow
    : public Indexable<T,MatrixRow<T,A,M,D>,Full,Abstract,VectorShape>
{
public:
    typedef Indexable<T,MatrixRow<T,A,M,D>,Full,Abstract,VectorShape> IndexableT;
    typedef Ref<T,IndexableT,VectorShape> RefT;

    MatrixRow(Indexable<T,A,M,D,MatrixShape>& m,index_t row)
        : itsMatrix(static_cast<A&>(m))
        , itsRow   (row)
    {};
    MatrixRow(const Indexable<T,A,M,D,MatrixShape>& m,index_t row)
        : itsMatrix(static_cast<A&>(const_cast<Indexable<T,A,M,D,MatrixShape>&>(m)))
        , itsRow   (row)
    {};
    MatrixRow(const MatrixRow& r)
        : itsMatrix(r.itsMatrix)
        , itsRow   (r.itsRow)
    {};
    ~MatrixRow() {};
    MatrixRow& operator=(const MatrixRow& r)
    {
        VectorAssign(*this,r);
        return *this;
    }
    MatrixRow& operator=(T t)
    {
        VectorAssign(*this,t);
        return *this;
    }
    template <class A1,Store M1, Data D1> MatrixRow& operator=(const Indexable<T,A1,M1,D1,VectorShape>& v)
    {
        VectorAssign(*this,v);
        return *this;
    }

    index_t   size  () const
    {
        return itsMatrix.GetLimits().Col.size();
    }
    VecLimits GetLimits() const
    {
        return itsMatrix.GetLimits().Col;
    };

    T  operator()(index_t col) const
    {
        return itsMatrix(itsRow,col);
    }
    T& operator()(index_t col)
    {
        return itsMatrix(itsRow,col);
    }

    MatrixRow& operator*=(const T& t)
    {
        for (index_t c:itsMatrix.cols())
            (*this)(c)*=t;
        return *this;
    }
    MatrixRow& operator/=(const T& t)
    {
        for (index_t c:itsMatrix.cols())
            (*this)(c)/=t;
        return *this;
    }

    class Subscriptor
    {
    public:
        Subscriptor(Indexable<T,MatrixRow<T,A,M,D>,Full,Abstract,VectorShape>& a)
            : itsRow(static_cast<MatrixRow<T,A,M,D>&>(a).itsRow)
            , itsS(static_cast<MatrixRow<T,A,M,D>&>(a).itsMatrix)
        {};

        T& operator()(index_t i)
        {
            return itsS(itsRow,i);    //Let matrix do bounds check.
        }

    private:
        const index_t itsRow;
        typename A::Subscriptor itsS;
    };

private:
    friend class Indexable<T,MatrixRow,Full,Abstract,VectorShape>;
    friend class Subscriptor;

    A&    itsMatrix;
    const index_t itsRow;
};



//---------------------------------------------------------------------
//
//  Matrix column proxy class.
//
template <class T, class A, Store M, Data D> class MatrixColumn
    : public Indexable<T,MatrixColumn<T,A,M,D>,Full,Abstract,VectorShape>
{
public:
    typedef Indexable<T,MatrixColumn<T,A,M,D>,Full,Abstract,VectorShape> IndexableT;
    typedef Ref<T,IndexableT,VectorShape> RefT;
    MatrixColumn(Indexable<T,A,M,D,MatrixShape>& m,index_t col)
        : itsMatrix(static_cast<A&>(m))
        , itsColumn(col)
    {};
    MatrixColumn(const Indexable<T,A,M,D,MatrixShape>& m,index_t col)
        : itsMatrix(static_cast<A&>(const_cast<Indexable<T,A,M,D,MatrixShape>&>(m)))
        , itsColumn(col)
    {};
    MatrixColumn(const MatrixColumn& r)
        : itsMatrix(r.itsMatrix)
        , itsColumn(r.itsColumn)
    {};
    ~MatrixColumn() {};
    MatrixColumn& operator=(const MatrixColumn& r)
    {
        VectorAssign(*this,r);
        return *this;
    }
    MatrixColumn& operator=(T t)
    {
        VectorAssign(*this,t);
        return *this;
    }
    template <class A1,Store M1, Data D1> MatrixColumn& operator=(const Indexable<T,A1,M1,D1,VectorShape>& v)
    {
        VectorAssign(*this,v);
        return *this;
    }

    index_t   size  () const
    {
        return itsMatrix.GetLimits().Row.size();
    }
    VecLimits GetLimits() const
    {
        return itsMatrix.GetLimits().Row;
    };

    T  operator()(index_t row) const
    {
        return itsMatrix(row,itsColumn);
    }
    T& operator()(index_t row)
    {
        return itsMatrix(row,itsColumn);
    }

    MatrixColumn& operator*=(const T& t)
    {
        for (index_t r:itsMatrix.rows())
            (*this)(r)*=t;
        return *this;
    }
    MatrixColumn& operator/=(const T& t)
    {
        for (index_t r:itsMatrix.rows())
            (*this)(r)/=t;
        return *this;
    }


    class Subscriptor
    {
    public:
        Subscriptor(Indexable<T,MatrixColumn<T,A,M,D>,Full,Abstract,VectorShape>& a)
            : itsColumn(static_cast<MatrixColumn<T,A,M,D>&>(a).itsColumn)
            , itsS(static_cast<MatrixColumn<T,A,M,D>&>(a).itsMatrix)
        {};

        T& operator()(index_t i)
        {
            return itsS(i,itsColumn);    //Let matrix do bounds check.
        }

    private:
        const index_t itsColumn;
        typename A::Subscriptor itsS;
    };

private:
    friend class Indexable<T,MatrixColumn,Full,Abstract,VectorShape>;
    friend class Subscriptor;

    A&    itsMatrix;
    const index_t itsColumn;
};


//---------------------------------------------------------------------
//
//  Matrix diagonal proxy class.
//
template <class T, class A, Store M, Data D> class MatrixDiagonal
    : public Indexable<T,MatrixDiagonal<T,A,M,D>,Full,Abstract,VectorShape>
{
public:
    typedef Indexable<T,MatrixDiagonal<T,A,M,D>,Full,Abstract,VectorShape> IndexableT;
     typedef Ref<T,IndexableT,VectorShape> RefT;

    MatrixDiagonal(Indexable<T,A,M,D,MatrixShape>& m)
        : itsMatrix(static_cast<A&>(m))
    {
//        assert(m.GetLimits().Row==m.GetLimits().Col);
    };
    MatrixDiagonal(const Indexable<T,A,M,D,MatrixShape>& m)
        : itsMatrix(static_cast<A&>(const_cast<Indexable<T,A,M,D,MatrixShape>&>(m)))
    {
//        assert(m.GetLimits().Row==m.GetLimits().Col);
    };
    MatrixDiagonal(const MatrixDiagonal& m)
        : itsMatrix(m.itsMatrix)
    {};
    ~MatrixDiagonal() {};
    MatrixDiagonal& operator=(const MatrixDiagonal& r)
    {
        VectorAssign(*this,r);
        return *this;
    }
    MatrixDiagonal& operator=(T t)
    {
        VectorAssign(*this,t);
        return *this;
    }
    template <class A1,Store M1, Data D1> MatrixDiagonal& operator=(const Indexable<T,A1,M1,D1,VectorShape>& v)
    {
        VectorAssign(*this,v);
        return *this;
    }

    index_t   size  () const
    {
        return std::min(itsMatrix.GetLimits().Row.size(),itsMatrix.GetLimits().Col.size());
    }
    VecLimits GetLimits() const
    {
        return itsMatrix.GetLimits().Row.size() >= itsMatrix.GetLimits().Col.size()
            ? itsMatrix.GetLimits().Col
            : itsMatrix.GetLimits().Row;
    };

    T  operator()(index_t i) const
    {
        return itsMatrix(i,i);
    }
    T& operator()(index_t i)
    {
        return itsMatrix(i,i);
    }

    class Subscriptor
    {
    public:
        Subscriptor(Indexable<T,MatrixDiagonal<T,A,M,D>,Full,Abstract,VectorShape>& a)
            : itsS(static_cast<MatrixDiagonal<T,A,M,D>&>(a).itsMatrix)
        {};

        T& operator()(index_t i)
        {
            return itsS(i,i);    //Let matrix do bounds check.
        }

    private:
        typename A::Subscriptor itsS;
    };

private:
    friend class Indexable<T,MatrixDiagonal,Full,Abstract,VectorShape>;
    friend class Subscriptor;

    A&    itsMatrix;
};

//----------------------------------------------------------------------------
//
//  Full matrix class with conventional direct subscripting, i.e. no list of
//  column pointers is maintained.
//
export template <class T> class Matrix
  : public MatrixBase
  , public ArrayIndexable<T,Matrix<T>,Full     ,MatrixShape>
  , public Indexable     <T,Matrix<T>,Full,Real,MatrixShape>
  , public TStreamableObject<Matrix<T> >
{
 public:
  typedef ArrayIndexable<T,Matrix<T>,Full     ,MatrixShape> IterableT;
  typedef      Indexable<T,Matrix<T>,Full,Real,MatrixShape> IndexableT;
  typedef Ref<T,IndexableT,MatrixShape> RefT;

  explicit Matrix(                );
  explicit Matrix( size_t  r, size_t c);
  explicit Matrix(index_t rl,index_t rh,index_t cl,index_t ch);
  explicit Matrix(const VecLimits& r,const VecLimits& c);
  explicit Matrix(const MatLimits&);
  Matrix(const Matrix& m);
  Matrix(Matrix& m) : Matrix<T>(const_cast<const Matrix&>(m)) {};
  template <class A>         Matrix(const ArrayIndexable<T,A,Full,MatrixShape>&);
  template <class A,Store M> Matrix(const Indexable<T,A,M,Abstract,MatrixShape>&);
  template <class A,Store M> Matrix(const Indexable<T,A,M,Real    ,MatrixShape>&);

  Matrix& operator=(const Matrix&);
  template <class A>         Matrix& operator=(const ArrayIndexable<T,A,Full,MatrixShape>&);
  template <class A,Store M> Matrix& operator=(const Indexable<T,A,M,Abstract,MatrixShape>&);
  template <class A,Store M> Matrix& operator=(const Indexable<T,A,M,Real    ,MatrixShape>&);

#ifdef OML_MOVE_OPS
  Matrix(Matrix&& m);
  Matrix& operator=(Matrix&&);
#endif

  MatLimits ReBase(index_t rlow,index_t clow);
  MatLimits ReBase(const MatLimits& lim);


  std::ostream& Write(std::ostream&) const;
  std::istream& Read (std::istream&)      ;
  friend std::ostream& operator<<(std::ostream& os,const Matrix& a)
  {
    return os << static_cast<const TStreamableObject<Matrix<T> >& >(a);
  }

  const T& ref(index_t i,index_t j) const;
  const T& operator()(index_t,index_t) const;
        T& operator()(index_t,index_t)      ;

  size_t    size  () const; //Required by iterable.
  MatLimits GetLimits() const;

  void SetLimits(const MatLimits&                           , bool preserve=false);
  void SetLimits(size_t  r , size_t c                       , bool preserve=false);
  void SetLimits(index_t rl,index_t rh,index_t cl,index_t ch, bool preserve=false);
  void SetLimits(const VecLimits& r,const VecLimits& c      , bool preserve=false);
  void RemoveRow   (index_t r);
  void RemoveColumn(index_t r);

  void ReIndexRows   (const std::vector<index_t>& index);
  void ReIndexColumns(const std::vector<index_t>& index);
  void SwapRows   (index_t i,index_t j);
  void SwapColumns(index_t i,index_t j);
  Matrix SubMatrix(const MatLimits& lim) const;

  typedef MatrixRow     <T,Matrix<T>,Full,Real> RowType;
  typedef MatrixColumn  <T,Matrix<T>,Full,Real> ColType;
  typedef MatrixDiagonal<T,Matrix<T>,Full,Real> DiagType;

  RowType  GetRow     (index_t row) {return RowType (*this,row);}
  ColType  GetColumn  (index_t col) {return ColType (*this,col);}
  DiagType GetDiagonal(           ) {return DiagType(*this    );}

  const RowType  GetRow     (index_t row) const {return RowType (*this,row);}
  const ColType  GetColumn  (index_t col) const {return ColType (*this,col);}
  const DiagType GetDiagonal(           ) const {return DiagType(*this    );}

  typedef typename IterableT::const_iterator  const_iterator ;
  typedef typename IterableT::iterator iterator;

  Matrix& operator+=(const ArrayIndexable<T,Matrix,Full,MatrixShape>& b) {return ArrayAdd(*this,b);}
  Matrix& operator-=(const ArrayIndexable<T,Matrix,Full,MatrixShape>& b) {return ArraySub(*this,b);}
  template <class B> Matrix& operator+=(const Indexable<T,B,Full,Abstract,MatrixShape>& b)
  {
      if (size()==0)
      {
        SetLimits(b.GetLimits(),false);
        Fill(*this,T(0));
      }
      return MatrixAdd(*this,b);
  }
  template <class B> Matrix& operator-=(const Indexable<T,B,Full,Abstract,MatrixShape>& b)
  {
      if (size()==0)
      {
        SetLimits(b.GetLimits(),false);
        Fill(*this,T(0));
      }
      return MatrixSub(*this,b);
  }


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
      , itsPtr(static_cast<Matrix&>(m).priv_begin())
      {};

    T& operator()(index_t i,index_t j) {CHECK(i,j);return itsPtr[itsLimits.Offset(i,j)];}

  private:
    MatLimits itsLimits;
    T*        itsPtr;
  };
#undef CHECK

 private:
  friend class Indexable<T,Matrix,Full,Real,MatrixShape>;
  friend class ArrayIndexable <T,Matrix,Full     ,MatrixShape>;
  friend class Subscriptor;

  const T* priv_begin() const {return &*itsData.begin();} //Required by ArrayIndexable.
        T* priv_begin()       {return &*itsData.begin();} //Required by ArrayIndexable.
  void  Check () const; //Check internal consistency between limits and cow.

#ifdef OML_USE_STDVEC
   std::vector<T> itsData;
#else
  cow_array<T> itsData;   //Copy-On-Write array for the data.
#endif
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


template <class T> inline const T& Matrix<T>::operator()(index_t i,index_t j) const
{
  CHECK(i,j);
  return itsData[GetLimits().Offset(i,j)];
}

template <class T> inline T& Matrix<T>::operator()(index_t i,index_t j)
{
  CHECK(i,j);
  return itsData[GetLimits().Offset(i,j)];
}

template <class T> inline  const T& Matrix<T>::ref(index_t i,index_t j) const
{
  CHECK(i,j);
  index_t index=GetLimits().Offset(i,j);
  return itsData[index];
}

#undef CHECK


#ifdef DEBUG
#define CHECK Check()
#else
#define CHECK
#endif

//-----------------------------------------------------------------------------
//
//  Constructors.
//
template <class T> Matrix<T>::Matrix()
  : itsData(0)
  {
    CHECK;
  }

template <class T> Matrix<T>::Matrix(size_t r, size_t c)
  : MatrixBase(MatLimits(r,c))
  , itsData   (GetLimits().size()   )
  {
    CHECK;
  }

template <class T> Matrix<T>::Matrix(index_t rl,index_t rh,index_t cl,index_t ch)
  : MatrixBase(MatLimits(rl,rh,cl,ch))
  , itsData   (GetLimits().size()   )
  {
    CHECK;
  }

template <class T> Matrix<T>::Matrix(const VecLimits& r,const VecLimits& c)
  : MatrixBase(MatLimits(r,c))
  , itsData   (GetLimits().size()   )
  {
    CHECK;
  }

template <class T> Matrix<T>::Matrix(const MatLimits& lim)
  : MatrixBase(lim          )
  , itsData   (lim.size())
  {
    CHECK;
  }

template <class T> Matrix<T>::Matrix(const Matrix<T>& m)
  : MatrixBase(m)
  , itsData   (m.itsData)
  {
//    std::cout << "Matrix<T> copy constructor" << std::endl;
    CHECK;
  }


template <class T> Matrix<T>& Matrix<T>::operator=(const Matrix<T>& m)
{
  MatrixBase::operator=(m);
  itsData=m.itsData;
//  std::cout << "Matrix<T> op=" << std::endl;
  return *this;
}


//-----------------------------------------------------------------------------
//
//  Internal consistency check.
//
template <class T> void Matrix<T>::Check() const
{
  assert(GetLimits().Check());
  assert(GetLimits().size()==itsData.size());
}


//-----------------------------------------------------------------------------
//
//  IO
//
template <class T> std::ostream& Matrix<T>::Write(std::ostream& os) const
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

template <class T> std::istream& Matrix<T>::Read(std::istream& is)
{
  assert(is);
  MatLimits lim;
  is >> lim;
  if (size()==0)
    SetLimits(lim);
   else
    assert(size()==0);

  ::Read(is,*this);
  CHECK;
  assert(is);
  return is;
}
//-----------------------------------------------------------------------------
//
//  Changing size and/or limits.
//
template <class T> void Matrix<T>::SetLimits(const MatLimits& theLimits, bool preserve)
{
#ifdef DEBUG
  theLimits.Check();
#endif
  if (GetLimits()!=theLimits)
  {
    if (preserve)
    {
        Matrix<T> dest(theLimits);

      index_t newrowlow=theLimits.Row.Low, newrowhigh=theLimits.Row.High;
      index_t newcollow=theLimits.Col.Low, newcolhigh=theLimits.Col.High;
      index_t rowlow =Max(newrowlow ,GetLimits().Row.Low );   //Limits of old and new
      index_t rowhigh=Min(newrowhigh,GetLimits().Row.High);   //data overlap.
      index_t collow =Max(newcollow ,GetLimits().Col.Low );   //Limits of old and new
      index_t colhigh=Min(newcolhigh,GetLimits().Col.High);   //data overlap.

      Subscriptor sdest (dest);         //Destination subscriptor.

      for (index_t i=rowlow;i<=rowhigh;i++)
        for (index_t j=collow;j<=colhigh;j++)
          sdest(i,j)=(*this)(i,j); //Transfer any overlaping data.

			(*this)=dest;
		}
     else
    {
			(*this)=Matrix<T>(theLimits);
    }
  }
}

//-----------------------------------------------------------------------------
//
//  Remove a row or a column
//
template <class T> void Matrix<T>::RemoveRow(index_t r)
{
    assert(r>=GetLimits().Row.Low);
    assert(r<=GetLimits().Row.High);
    MatLimits lim(GetLimits().Row.Low,GetLimits().Row.High-1,GetLimits().Col.Low,GetLimits().Col.High);

    Matrix<T> dest(lim);
    Subscriptor sdest (dest);         //Destination subscriptor.
    for (index_t j:dest.cols())
    {
        for (index_t i=lim.Row.Low; i<r; i++)
            sdest(i,j)=(*this)(i,j); //Copy low rows
        for (index_t i=r; i<=lim.Row.High; i++)
            sdest(i,j)=(*this)(i+1,j); //Shift high rows down by one
    }

    (*this)=dest;
}
template <class T> void Matrix<T>::RemoveColumn(index_t c)
{
    assert(c>=GetLimits().Col.Low);
    assert(c<=GetLimits().Col.High);
    MatLimits lim(GetLimits().Row.Low,GetLimits().Row.High,GetLimits().Col.Low,GetLimits().Col.High-1);

    Matrix<T> dest(lim);
    Subscriptor sdest (dest);         //Destination subscriptor.
    for (index_t i:dest.rows())
    {
        for (index_t j=lim.Col.Low; j<c; j++)
            sdest(i,j)=(*this)(i,j); //Copy low cols
        for (index_t j=c; j<=lim.Col.High; j++)
            sdest(i,j)=(*this)(i,j+1); //Shift high cols down by one
    }

    (*this)=dest;
}

template <class T> void Matrix<T>::ReIndexRows(const std::vector<index_t>& index)
{
  assert(GetLimits().GetNumRows()==index.size());

  typename std::vector<index_t>::const_iterator i=index.begin();
  Matrix<T> dest(GetLimits());
  Subscriptor sdest(dest);

  for (int row=GetLimits().Row.Low;row<=GetLimits().Row.High;row++,i++)
  for (int col=GetLimits().Col.Low;col<=GetLimits().Col.High;col++)
  sdest(row,col)=(*this)(*i+GetLimits().Row.Low,col);

  *this=dest;
}

template <class T> void Matrix<T>::ReIndexColumns(const std::vector<index_t>& index)
{
  assert(GetLimits().GetNumCols()==index.size());

  typename std::vector<index_t>::const_iterator i=index.begin();
  Matrix<T> dest(GetLimits());
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

template <class T> void Matrix<T>::SwapRows(index_t i,index_t j)
{
  Subscriptor s(*this);
  for (index_t c : this->cols())
  {
     T temp=s(i,c);
     s(i,c)=s(j,c);
     s(j,c)=temp;
  }
}

template <class T> void Matrix<T>::SwapColumns(index_t i,index_t j)
{
  Subscriptor s(*this);
  for (index_t r : this->rows())
  {
     T temp=s(r,i);
     s(r,i)=s(r,j);
     s(r,j)=temp;
  }
}

template <class T> Matrix<T> Matrix<T>::SubMatrix(const MatLimits& lim) const
{
    Matrix<T> dest(lim);
	assert(dest.GetLimits().Row.Low >=GetLimits().Row.Low );
	assert(dest.GetLimits().Col.Low >=GetLimits().Col.Low );
	assert(dest.GetLimits().Row.High<=GetLimits().Row.High);
	assert(dest.GetLimits().Col.High<=GetLimits().Col.High);
	Subscriptor s(dest);
	index_t rh=dest.GetLimits().Row.High;
	index_t ch=dest.GetLimits().Col.High;
	for (index_t i=dest.GetLimits().Row.Low;i<=rh;i++)
		for (index_t j=dest.GetLimits().Col.Low;j<=ch;j++)
			s(i,j)=(*this)(i,j);
    return dest;
}

#undef CHECK


template <class T> template <class A> inline
Matrix<T>::Matrix(const ArrayIndexable<T,A,Full,MatrixShape>& m)
  : MatrixBase(m.GetLimits        ())
  , itsData   (GetLimits().size())
  {
    ArrayAssign(*this,m); //Use op[].
  }

template <class T> template <class A,Store M> inline
Matrix<T>::Matrix(const Indexable<T,A,M,Abstract,MatrixShape>& m)
  : MatrixBase(m.GetLimits        ())
  , itsData   (GetLimits().size())
  {
    //m.GetLimits();
    MatrixAssign(*this,m); //Use op(i,j).
  }
template <class T> template <class A,Store M> inline
Matrix<T>::Matrix(const Indexable<T,A,M,Real,MatrixShape>& m)
  : MatrixBase(m.GetLimits        ())
  , itsData   (GetLimits().size())
  {
    //m.GetLimits();
    MatrixAssign(*this,m); //Use op(i,j).
  }


template <class T> template <class A> inline
Matrix<T>& Matrix<T>::operator=(const ArrayIndexable<T,A,Full,MatrixShape>& m)
{
  SetLimits(m.GetLimits());
  ArrayAssign(*this,m); //Use op[].
  return *this;
}

template <class T> template <class A,Store M> inline
Matrix<T>& Matrix<T>::operator=(const Indexable<T,A,M,Abstract,MatrixShape>& m)
{
  SetLimits(m.GetLimits());
  MatrixAssign(*this,m); //Use op(,).
  return *this;
}

template <class T> template <class A,Store M> inline
Matrix<T>& Matrix<T>::operator=(const Indexable<T,A,M,Real,MatrixShape>& m)
{
  SetLimits(m.GetLimits());
  MatrixAssign(*this,m); //Use op(,).
  return *this;
}

#ifdef OML_MOVE_OPS

template <class T> inline Matrix<T>::Matrix(Matrix<T>&& m)
  : MatrixBase(m)
  , itsData   (std::move(m.itsData))
  {
//    std::cout << "Matrix<T> move constructor m.itsData.size()=" << m.itsData.size() << std::endl;
  }

template <class T> inline Matrix<T>& Matrix<T>::operator=(Matrix<T>&& m)
{
  MatrixBase::operator=(m);
  itsData=std::move(m.itsData);
//  std::cout << "Matrix<T> move op=" << std::endl;
  return *this;
}
#endif

template <class T> inline MatLimits Matrix<T>::ReBase(index_t rlow,index_t clow)
{
    MatLimits oldLimits=GetLimits();
    MatrixBase::ReBase(rlow,clow);
    return oldLimits;
}

template <class T> inline MatLimits Matrix<T>::ReBase(const MatLimits& lim)
{
    MatLimits oldLimits=GetLimits();
    MatrixBase::ReBase(lim);
    return oldLimits;
}

template <class T, class B,Store M, Data D> Matrix<T>& operator*=(Matrix<T>& a,const Indexable<T,B,M,D,MatrixShape>& b)
{
    assert(a.GetLimits().Col==b.GetLimits().Row);
    Matrix<T> temp=a*b;
    a.SetLimits(temp.GetLimits());
    a=std::move(temp);
    return a;
}

template <class T> inline size_t  Matrix<T>::size() const
{
  return GetLimits().size();
}




template <class T> inline void Matrix<T>::SetLimits(size_t r, size_t c , bool preserve)
{
  SetLimits(MatLimits(r,c),preserve);
}

template <class T> inline void Matrix<T>::SetLimits(index_t rl,index_t rh,index_t cl,index_t ch, bool preserve)
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

// template class Matrix<double>;
//------------------------------------------------------------------
//
//  Matrix transpose proxy.
//
template <class T, class Mat> class MatrixTranspose
: public Indexable<T,MatrixTranspose<T,Mat>,Full,Abstract,MatrixShape>
{
 public:
  typedef Indexable<T,MatrixTranspose<T,Mat>,Full,Abstract,MatrixShape> IndexableT;
  typedef Ref<T,IndexableT,MatrixShape> RefT;

  MatrixTranspose(Mat m) : itsMatrix(m) {};
  MatrixTranspose(const MatrixTranspose<T,Mat>& mt) : itsMatrix(mt.itsMatrix) {};

  T operator()(index_t i, index_t j) const {return itsMatrix(j,i);}
  size_t  size() const {return GetLimits().size();}
  MatLimits GetLimits() const {return MatLimits(itsMatrix.GetLimits().Col,itsMatrix.GetLimits().Row);}

 private:
  Mat itsMatrix;
};

//-------------------------------------------------------------------------
//
//  Transpose operators, returns a proxy.
//
export template <class T, class A, Store M, Data D> inline
MatrixTranspose<T,Ref<T,Indexable<T,A,M,D,MatrixShape>,MatrixShape> >
Transpose(const Indexable<T,A,M,D,MatrixShape>& m)
{
	typedef Ref<T,Indexable<T,A,M,D,MatrixShape>,MatrixShape> ref;
  return MatrixTranspose<T,ref>(ref(m));
}

//! Overloaded operator version of Transpose, returns a proxy.
export template <class T, class A, Store M, Data D> inline
MatrixTranspose<T,Ref<T,Indexable<T,A,M,D,MatrixShape>,MatrixShape> >
operator~(const Indexable<T,A,M,D,MatrixShape>& m)
{
  typedef Ref<T,Indexable<T,A,M,D,MatrixShape>,MatrixShape> ref;
  return MatrixTranspose<T,ref>(ref(m));
}


//--------------------------------------------------------------------------
//
//  Matrix algebra functions  M*M, M*V V*M
//
//  All of these operators return a proxy for the operation.  So for example a transpose
//  proxy just stores a Matrix reference and op(i,j) just return ref(j,i). The multiplication
//  proxies are little more complicated, but essentially just return a Row(i)*Col(j) dot
//  product.



//----------------------------------------------------------------------
//
//  Multiplication, Matrix * Matrix Proxy.
//
template <class T, class A, class B> class MatrixMMOp
: public Indexable<T,MatrixMMOp<T,A,B>,Full,Abstract,MatrixShape>
{
 public:
  typedef Indexable<T,MatrixMMOp<T,A,B>,Full,Abstract,MatrixShape> IndexableT;
  typedef Ref<T,IndexableT,MatrixShape> RefT;

  MatrixMMOp(const A& a, const B& b)
    : itsA(a)
    , itsB(b)
  {
    //std::cout << "Alim="<< itsA.GetLimits() << " Blim=" << itsB.GetLimits()<< std::endl;
    assert(itsA.GetLimits().Col==itsB.GetLimits().Row);
  };
  MatrixMMOp(const MatrixMMOp& m)
    : itsA(m.itsA)
    , itsB(m.itsB)
  {}
  T operator()(index_t i, index_t j) const
  {
    T ret(0);
    for (index_t k:itsA.cols()) ret+=itsA(i,k)*itsB(k,j);
    return ret;
  }
  size_t  size() const {return GetLimits().size();}
  MatLimits GetLimits() const {return MatLimits(itsA.GetLimits().Row,itsB.GetLimits().Col);}

 private:
  const A itsA;
  const B itsB;
};

//----------------------------------------------------------------------
//
//  Multiplication, Matrix * Vector Proxy.
//
template <class T, class A, class V> class MatrixMVOp
: public Indexable<T,MatrixMVOp<T,A,V>,Full,Abstract,VectorShape>
{
 public:
  typedef Indexable<T,MatrixMVOp<T,A,V>,Full,Abstract,VectorShape> IndexableT;
  typedef Ref<T,IndexableT,VectorShape> RefT;

  MatrixMVOp(const A& a, const V& v)
    : itsA(a)
    , itsV(v)
  {
    assert(itsA.GetLimits().Col==itsV.GetLimits());
  };
  MatrixMVOp(const MatrixMVOp& m)
    : itsA(m.itsA)
    , itsV(m.itsV)
    {};
  T operator()(index_t i) const
  {
    T ret(0);
    for (index_t k:itsA.cols()) ret+=itsA(i,k)*itsV(k);
    return ret;
  }
  VecLimits GetLimits() const {return itsA.GetLimits().Row;}
  size_t    size  () const {return GetLimits().size();}

 private:
  const A itsA;
  const V itsV;
};

//----------------------------------------------------------------------
//
//  Multiplication, Vector * Matrix Proxy.
//
template <class T, class V, class B> class MatrixVMOp
: public Indexable<T,MatrixVMOp<T,V,B>,Full,Abstract,VectorShape>
{
 public:
  typedef Indexable<T,MatrixVMOp<T,V,B>,Full,Abstract,VectorShape> IndexableT;
  typedef Ref<T,IndexableT,VectorShape> RefT;

  MatrixVMOp(const V& v, const B& b)
    : itsV(v)
    , itsB(b)
  {
    assert(itsV.GetLimits()==itsB.GetLimits().Row);
  };
  MatrixVMOp(const MatrixVMOp& m)
    : itsV(m.itsV)
    , itsB(m.itsB)
    {};
  T operator()(index_t i) const
  {
    T ret(0);
    for (index_t k: itsB.rows()) ret+=itsV(k)*itsB(k,i);
    return ret;
  }
  VecLimits GetLimits() const {return itsB.GetLimits().Col;}
  size_t    size  () const {return GetLimits().size();}

 private:
  const V itsV;
  const B itsB;
};

//---------------------------------------------------------------------
//
//  Matrix * Matrix, returns a matrix instead of a proxy.
//
export template <class TA, class TB> inline
auto operator*(const Matrix<TA>& a,const Matrix<TB>& b)
{
  typedef typename ReturnType<TA,TB>::RetType TR;
  return Matrix<TR>(MatrixMMOp<TR,typename Matrix<TA>::RefT,typename Matrix<TB>::RefT>(a,b));
}
//---------------------------------------------------------------------
//
//  Matrix * Indexable<MatrixShape>, returns a matrix instead of a proxy.
//
export template <class TA, class TB, class B, Data DB> inline
auto operator*(const Matrix<TA>& a,const Indexable<TB,B,Full,DB,MatrixShape>& b)
{
  typedef typename ReturnType<TA,TB>::RetType TR;
  return Matrix<TR>(MatrixMMOp<TR,typename Matrix<TA>::RefT,typename B::RefT>(a,b));
}

//---------------------------------------------------------------------
//
//  Indexable<MatrixShape> * Matrix, returns a matrix.
//
export template <class TA, class TB, class A, Data DA> inline
auto operator*(const Indexable<TA,A,Full,DA,MatrixShape>& a,const Matrix<TB>& b)
{
  typedef typename ReturnType<TA,TB>::RetType TR;
  return Matrix<TR>(MatrixMMOp<TR,typename A::RefT,typename Matrix<TB>::RefT>(a,b));
}

//---------------------------------------------------------------------
//
//  Indexable<MatrixShape> * Indexable<MatrixShape>, returns a proxy.
//
export template <class TA, class TB, class A, class B, Data DA, Data DB> inline
auto operator*(const Indexable<TA,A,Full,DA,MatrixShape>& a,const Indexable<TB,B,Full,DB,MatrixShape>& b)
{
  typedef typename ReturnType<TA,TB>::RetType TR;
  return MatrixMMOp<TR,typename A::RefT,typename B::RefT>(a,b);
}


//--------------------------------------------------------------------------------
//
//  Matrix * Vector, returns a proxy.
//
export template <class TA, class TB, class A, class B, Store MA, Store MB, Data DA, Data DB> inline
auto operator*(const Indexable<TA,A,MA,DA,MatrixShape>& a,const Indexable<TB,B,MB,DB,VectorShape>& b)
{
  typedef typename ReturnType<TA,TB>::RetType TR;
  return MatrixMVOp<TR,typename A::RefT,typename B::RefT>(a,b);
}

//------------------------------------------------------------
//
//  Vector * Matrix, returns a proxy.
//
export template <class TA, class TB, class A, class B, Store MA, Store MB, Data DA, Data DB> inline
auto operator*(const Indexable<TA,A,MA,DA,VectorShape>& a,const Indexable<TB,B,MB,DB,MatrixShape>& b)
{
  typedef typename ReturnType<TA,TB>::RetType TR;
  return MatrixVMOp<TR,typename A::RefT,typename B::RefT>(a,b);
}

export template <class T> Matrix<T> TensorProduct(const Matrix<T>& a, const Matrix<T>& b)
{
    Matrix<T> r(a.GetLimits()*b.GetLimits());
    //
    //  These loops need to run in the same order as SiteOperatorImp::Product
    //
    int i=r.GetRowLimits().Low;
    for (index_t ib:b.rows())
    for (index_t ia:a.rows())
    {
        int j=r.GetColLimits().Low;
        for (index_t jb:b.cols())
        for (index_t ja:a.cols())
        {
            r(i,j)=a(ia,ja)*b(ib,jb);
            j++;
        }
        i++;
    }
    return r;
}

//--------------------------------------------------------------
//
//  Other assorted matrix functions.
//

// Create a unit of Kronecker Delta matrix.
export template <class T, class A, Store M> inline
void Unit(Indexable<T,A,M,Real,MatrixShape>& m)
{
  A& a(static_cast<A&>(m));
  Fill(a,T(0.0));
  a.GetDiagonal()=T(1);
}

export template <class T, class A, Store M, Data D> inline
A operator~(const Indexable<std::complex<T>,A,M,D,MatrixShape>& m)
{
  return conj(Transpose(m));
}

// Check if matrix is symmetric. This allows for no roundoff errors.
export template <class T, class A,Data D> inline
bool IsSymmetric(const Indexable<T,A,Full,D,MatrixShape>& m)
{
  return m==Transpose(m);
}

// Check if matrix is symmetric. This allows for some tolerance.
export template <class T, class A,Data D> inline
bool IsSymmetric(const Indexable<T,A,Full,D,MatrixShape>& m,double eps)
{
  return Max(fabs(m-Transpose(m)))<=eps;
}

// Check if matrix is anit-symmetric.
export template <class T, class A, Data D> inline
bool IsAntiSymmetric(const Indexable<T,A,Full,D,MatrixShape>& m)
{
  return m==-Transpose(m);
}

// Check if complex matrix  is Hermitian.
export template <class T, class A, Data D> inline
bool IsHermitian(const Indexable<std::complex<T>,A,Full,D,MatrixShape>& m)
{
  return m==conj(Transpose(m));
}

// Check if complex matrix  is Hermitian to within a specified precision eps.
export template <class T, class A, Data D> inline
bool IsHermitian(const Indexable<std::complex<T>,A,Full,D,MatrixShape>& m, double eps)
{
  return Max(fabs(m-conj(Transpose(m))))<=eps;
}

// Check if complex matrix  is normal to within a specified precision eps.
export template <class T, class A, Data D> inline
bool IsNormal(const Indexable<T,A,Full,D,MatrixShape>& m, double eps)
{
  Matrix<T> mdagger=~m;
  return Max(fabs(m*mdagger-mdagger*m))<=eps;
}


// Dummy check for real matrices.
export template <class T, class A, Data D> inline
bool IsHermitian(const Indexable<T,A,Full,D,MatrixShape>& m, double eps)
{
  return IsSymmetric(m,eps);
}

// Force a Matrix to be symmetric by averaging off diagonal elements.
export template <class T, class A>
double MakeSymmetric(Indexable<T,A,Full,Real,MatrixShape>& m)
{
  A temp=Transpose(m);
  double del=Max(fabs(m-temp));
  static_cast<A&>(m)=(m+temp)/T(2);
  return del;
}

export template <class T, class A, Data D> inline
double FrobeniusNorm(const Indexable<T,A,Full,D,MatrixShape>& m)
{
    double fnorm=0.0;
    for (index_t i: m.rows())
        for (index_t j: m.cols())
            fnorm+=real(m(i,j)*conj(m(i,j)));
    return sqrt(fnorm);
}

export template <class T, class A, Data D> inline
bool IsUnit(const Indexable<T,A,Full,D,MatrixShape>& m)
{
    Matrix<T> U(m.GetLimits());
    Unit(U);
    return Max(fabs(m-U))==0.0;
}

export template <class T, class A, Data D> inline
bool IsUnit(const Indexable<T,A,Full,D,MatrixShape>& m,double eps)
{
    Matrix<T> U(m.GetLimits());
    Unit(U);
    return Max(fabs(m-U))<=eps;
}

export template <class T, class A, Data D> inline
bool IsDiagonal(const Indexable<T,A,Full,D,MatrixShape>& m,double eps)
{
    double off(0.0);
    for (index_t i:m.rows())
        for (index_t j:m.cols())
            if (i!=j) off=Max(off,fabs(m(i,j)));
    return off<=eps;
}

export template <class T, class A, Data D> inline
bool IsDiagonal(const Indexable<T,A,Full,D,MatrixShape>& m)
{
    return IsDiagonal(m,0.0);
}

export template <class T, class A, Data D> inline
bool IsLowerTriangular(const Indexable<T,A,Full,D,MatrixShape>& m)
{
    bool ret=true;
    MatLimits l=m.GetLimits();
    size_t Nr=l.GetNumRows(),Nc=l.GetNumCols();
    int delta=Nc>Nr ? Nc-Nr : 0;
    for (index_t i: m.rows())
        for (index_t j:m.cols(i+delta+1))
        {
            ret = ret && (m(i,j)==0.0);
            if (!ret) break;
        }
    return ret;
}

export template <class T, class A, Data D> inline
bool IsLowerTriangular(const Indexable<T,A,Full,D,MatrixShape>& m,double eps)
{
    bool ret=true;
    MatLimits l=m.GetLimits();
    size_t Nr=l.GetNumRows(),Nc=l.GetNumCols();
    int delta=Nc>Nr ? Nc-Nr : 0;
    for (index_t i: m.rows())
        for (index_t j:m.cols(i+delta+1))
        {
            ret = ret && (fabs(m(i,j))<=eps);
            if (!ret) break;
        }
    return ret;
}

export template <class T, class A, Data D> inline
bool IsUpperTriangular(const Indexable<T,A,Full,D,MatrixShape>& m)
{
    bool ret=true;
    MatLimits l=m.GetLimits();
    size_t Nr=l.GetNumRows(),Nc=l.GetNumCols();
    int delta=Nr>Nc ? Nr-Nc : 0;
    if (m.GetLimits().GetNumRows()!=0 && m.GetLimits().GetNumCols()!=0)
        for (index_t j: m.cols())
            for (index_t i:m.rows(j+delta+1))
            {
                ret = ret && (m(i,j)==0.0);
                if (!ret) break;
            }
    return ret;
}

export template <class T, class A, Data D> inline
bool IsUpperTriangular(const Indexable<T,A,Full,D,MatrixShape>& m,double eps)
{
    bool ret=true;
    MatLimits l=m.GetLimits();
    size_t Nr=l.GetNumRows(),Nc=l.GetNumCols();
    int delta=Nr>Nc ? Nr-Nc : 0;
    if (m.GetLimits().GetNumRows()!=0 && m.GetLimits().GetNumCols()!=0)
        for (index_t j: m.cols())
            for (index_t i:m.rows(j+delta+1))
            {
                ret = ret && (fabs(m(i,j))<=eps);
                if (!ret) break;
            }
    return ret;
}

export template <class T, class A, Data D> bool IsTriangular(Store ul,const Indexable<T,A,Full,D,MatrixShape>& m,double eps)
{
    bool ret=true;
    if (m.GetLimits().GetNumCols()>1 && m.GetLimits().GetNumRows()>1)
    {
        switch (ul)
        {
        case Upper:
            ret=IsUpperTriangular(m,eps);
            break;
        case Lower:
            ret=IsLowerTriangular(m,eps);
            break;
        default:
            ret=false;
        }
    }
    return ret;
}
//----------------------------------------------------------------------
//
//  Vector outer product Proxy.
//
template <class T,class V1,class V2> class MatrixOuterProduct
: public Indexable<T,MatrixOuterProduct<T,V1,V2>,Full,Abstract,MatrixShape>
{
 public:
  typedef Indexable<T,MatrixOuterProduct<T,V1,V2>,Full,Abstract,MatrixShape> IndexableT;
  typedef Ref<T,IndexableT,MatrixShape> RefT;

  MatrixOuterProduct(const V1& a, const V2& b) : itsA(a), itsB(b) {};
  MatrixOuterProduct(const MatrixOuterProduct& op) : itsA(op.itsA), itsB(op.itsB) {};

  T operator()(index_t i, index_t j) const {return itsA(i)*itsB(j);}
  MatLimits GetLimits() const {return MatLimits(itsA.GetLimits(),itsB.GetLimits());}
  size_t  size() const {return MatLimits().size();}

 private:
  const V1  itsA;
  const V2  itsB;
};

template <class T,class V> class SymMatrixOuterProduct
: public Indexable<T,SymMatrixOuterProduct<T,V>,Symmetric,Abstract,MatrixShape>
{
 public:
  typedef Indexable<T,SymMatrixOuterProduct<T,V>,Symmetric,Abstract,MatrixShape> IndexableT;
  typedef Ref<T,IndexableT,MatrixShape> RefT;

  SymMatrixOuterProduct(const V& v) : itsV(v) {};
  SymMatrixOuterProduct(const SymMatrixOuterProduct& op) : itsV(op.itsV) {};

  T operator()(index_t i, index_t j) const {return itsV(i)*itsV(j);}
  MatLimits GetLimits() const {return MatLimits(itsV.GetLimits(),itsV.GetLimits());}
  size_t  size() const {return MatLimits().size();}

 private:
  const V itsV;
};

//--------------------------------------------------------------------------
//
//  Outer product, returns a proxy.
//
export template <class T, class A, class B, Store MA, Store MB, Data DA, Data DB> inline
MatrixOuterProduct<T,Ref<T,Indexable<T,A,MA,DA,VectorShape>,VectorShape>,Ref<T,Indexable<T,B,MB,DB,VectorShape>,VectorShape> >
OuterProduct(const Indexable<T,A,MA,DA,VectorShape>& v1,const Indexable<T,B,MB,DB,VectorShape>& v2)
{
  typedef typename A::RefT refa;
  typedef typename B::RefT refb;
  return MatrixOuterProduct<T,refa,refb>(refa(v1),refb(v2));
}

export template <class T, class A, Store M, Data D> inline
SymMatrixOuterProduct<T,Ref<T,Indexable<T,A,M,D,VectorShape>,VectorShape> >
OuterProduct(const Indexable<T,A,M,D,VectorShape>& v)
{
  typedef Ref<T,Indexable<T,A,M,D,VectorShape>,VectorShape> ref;
  return SymMatrixOuterProduct<T,ref>(ref(v));
}



