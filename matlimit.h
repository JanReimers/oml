// File: matlimit.H  Vector and matrix limits structures, header file.
#ifndef _matlimit_H_
#define _matlimit_H_

// Copyright (1994-2005), Jan N. Reimers

#include "oml/veclimit.h"

//---------------------------------------------------------------------------
/*! \class MatLimits matlimit.h oml/matlimit.h
  \brief Encapsulate all effects of matrix subscript bounds in this class.

  Structure for matrix dimension limits.  All limits are checked in the
  VecLimits constructors.

  \nosubgrouping
*/
class MatLimits
{
 public:
  /*! \name Constructors*/
  //@{
  MatLimits(                                 );//!<Limits for a null Matrix.
  MatLimits(subsc_t,subsc_t,subsc_t,subsc_t  );//!<Construct from lower and upper bounds.
  MatLimits(index_t,index_t                  );//!<Construct from row and column size, use default lower bound.
  MatLimits(const VecLimits&,const VecLimits&);//!<Construct from lower and upper bounds.
  //@}
 ~MatLimits();

  static index_t size(index_t,index_t                  );
  static index_t size(subsc_t,subsc_t,subsc_t,subsc_t  );
  static index_t size(const VecLimits&,const VecLimits&);

  //! Returns number of rows.
  index_t GetNumRows        (       ) const;
  //! Returns number of columns.
  index_t GetNumCols        (       ) const;
  //! Returns number of elements.
  index_t size           (       ) const;

  index_t Offset    (subsc_t,subsc_t) const;
  bool    Check     (               ) const;
  bool    CheckIndex(subsc_t,subsc_t) const;

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

inline MatLimits::MatLimits(subsc_t rowLow, subsc_t rowHigh,subsc_t colLow, subsc_t colHigh)
  : Row(rowLow,rowHigh)
  , Col(colLow,colHigh)
  {}

inline MatLimits::MatLimits(index_t rowSize, index_t colSize)
  : Row(rowSize)
  , Col(colSize)
  {}

inline MatLimits::MatLimits(const VecLimits& rowLimits, const VecLimits& colLimits)
  : Row(rowLimits)
  , Col(colLimits)
  {}

inline MatLimits::~MatLimits() {}

//------------------------------------------------------------------------
//
//  Get size static functions.
//
inline index_t MatLimits::size(index_t numRows, index_t numCols)
{
  return numRows*numCols;
}
inline index_t MatLimits::size(const VecLimits& row, const VecLimits& col)
{
  return row.size()*col.size();
}

inline index_t MatLimits::size(subsc_t rowLow, subsc_t rowHigh,subsc_t colLow, subsc_t colHigh)
{
  return size(VecLimits(rowLow,rowHigh),VecLimits(colLow,colHigh) );
}

inline index_t MatLimits::GetNumRows() const
{
  return Row.size();
}

inline index_t MatLimits::GetNumCols() const
{
  return Col.size();
}

inline index_t MatLimits::size() const
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

inline index_t MatLimits::Offset(subsc_t i,subsc_t j) const
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

#endif //_matlimits_H_
