module;
#include <cstddef>

export module oml.Matrixbase;
import oml.VecLimits;
import oml.MatLimits;

//----------------------------------------------------------------------------
/*! \class MatrixBase matrixbase.h oml/matrixbase.h
  \brief Provids all sorts of information on the row and column limits.
*/
export class MatrixBase
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
