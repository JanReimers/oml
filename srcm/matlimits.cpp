module;
#include <iostream>
#include <cassert>


export module oml.MatLimits;
import oml.VecLimits;
import oml.StreamableObject;


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
// IO
//
std::ostream& MatLimits::Write(std::ostream& os) const
{
  assert(os);
  if(StreamableObject::Binary()) os << Row << Col;
  if(StreamableObject::Ascii ()) os << Row << " " << Col << " ";
  if(StreamableObject::Pretty()) os << Row << "," << Col << " ";
  assert(os);
  return os;
}
std::istream& MatLimits::Read(std::istream& is)
{
  assert(is);
  if(StreamableObject::Binary()) is >> Row >> Col;
  if(StreamableObject::Ascii ()) is >> Row >> Col;
  assert(is);
  return is;
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

export inline std::ostream& operator<<(std::ostream& os,const MatLimits& lim)
{
  return lim.Write(os);
}

export inline std::istream& operator>>(std::istream& is, MatLimits& lim)
{
  return lim.Read (is);
}

// Use this to get limits of a tensor product
export inline MatLimits operator*(const MatLimits& a, const MatLimits& b)
{
    return MatLimits(a.Row*b.Row,a.Col*b.Col);
}
