// File: MatrixBase.H  Abstract matrix with limits only.
#ifndef _MatrixBase_H_
#define _MatrixBase_H_

// Copyright (1994-2003), Jan N. Reimers

#include "oml/matlimit.h"

//----------------------------------------------------------------------------
/*! \class MatrixBase matrixbase.h oml/matrixbase.h
  \brief Provids all sorts of information on the row and column limits.
*/
class MatrixBase
{
 public:
  MatrixBase(                    );
  MatrixBase(const MatLimits& lim);

  /*! \name Limits info
  */
  //@{
  //! Extract info on Matrix size/shape.
  MatLimits GetLimits   () const;
  VecLimits GetRowLimits() const;
  VecLimits GetColLimits() const;
  index_t   GetNumRows  () const;
  index_t   GetNumCols  () const;
  index_t   size        () const;
  subsc_t   GetRowLow   () const;
  subsc_t   GetRowHigh  () const;
  subsc_t   GetColLow   () const;
  subsc_t   GetColHigh  () const;
  //@}
  void SetLimits(const MatLimits&, bool preserve=false);

 private:
  MatLimits    itsLimits; //Manages the upper and lower std::vector limits.
};

inline MatrixBase::MatrixBase(                    ) : itsLimits(   ) {}
inline MatrixBase::MatrixBase(const MatLimits& lim) : itsLimits(lim) {}

inline void MatrixBase::SetLimits(const MatLimits& lim,bool /*preserve*/)
{
  itsLimits=lim;
}

inline  MatLimits MatrixBase::GetLimits() const
{
  return itsLimits;
}

inline index_t MatrixBase::GetNumRows() const
{
  return itsLimits.Row.size();
}

inline index_t MatrixBase::GetNumCols() const
{
  return itsLimits.Col.size();
}

inline index_t MatrixBase::size() const
{
  return itsLimits.size();
}

inline VecLimits MatrixBase::GetRowLimits() const
{
  return itsLimits.Row;
}

inline VecLimits MatrixBase::GetColLimits() const
{
  return itsLimits.Col;
}

inline subsc_t MatrixBase::GetRowLow() const
{
  return itsLimits.Row.Low;
}

inline subsc_t MatrixBase::GetRowHigh() const
{
  return itsLimits.Row.High;
}

inline subsc_t MatrixBase::GetColLow() const
{
  return itsLimits.Col.Low;
}

inline subsc_t MatrixBase::GetColHigh() const
{
  return itsLimits.Col.High;
}


#endif //_MatrixBase_H_
