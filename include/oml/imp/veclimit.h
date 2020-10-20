// File: veclimit.H  Vector and matrix limits data structures.
#ifndef _veclimit_H_
#define _veclimit_H_

// Copyright (1994-2005), Jan N. Reimers

#include "oml/imp/index_t.h"
#include <iosfwd>

//-----------------------------------------------------------------------------
//
//  Default upper and lower bounds.  This where the Fortan lower bound
//  of 1 gets decided.  Just make LOW=0 if want the C default indexing.
//
const index_t LOW=1;
const index_t HIGH=LOW-1;
//-----------------------------------------------------------------------------
/*! \class VecLimits veclimit.h oml/veclimit.h
  \brief Encapsulate all effects of non-zero lower subscript bounds in this class.

  This is just a data structure so the data is public and there is no
  attempt at any encapsulation stuff here, don't need it.
  \nosubgrouping
*/
class VecLimits
{
 public:
  /*! \name Constructors*/
  //@{
  VecLimits(               ); //!<Limits for a null Vector.
  VecLimits(index_t        ); //!<Construct from size, use default lower bound.
  VecLimits(index_t,index_t); //!<Construct from lower and upper bounds.
  //@}
 ~VecLimits();

  static index_t size(index_t,index_t)      ;
  //! Returns number of elements.
         index_t size(               ) const;

  index_t Offset         (index_t) const;
  bool    Check          (       ) const;
  bool    CheckIndex     (index_t) const;
  //! Comparison.
  bool operator==(const VecLimits&) const;
  //! Comparison.
  bool operator!=(const VecLimits&) const;

  /*! \name IO
    See StreamableObject for details on output formats. You must include std::vector_io.h to get
    definitions.  op>> and op<< are also defined.
  */
  //@{
  std::ostream& Write(std::ostream&) const; //!<Write to stream.
  std::istream& Read (std::istream&)      ; //!<Read from stream.
  //@}

  index_t Low;   //!<lower std::vector index limit.
  index_t High;  //!<upper std::vector index limit.
};

#ifdef DEBUG
#define CHECK Check()
#else
#define CHECK
#endif

//-----------------------------------------------------------------------------
//
//  Constructors.
//
inline VecLimits::VecLimits() :                 //Limits for a null std::vector.
  Low(LOW),
  High(HIGH)
  {
    CHECK;
  }

inline VecLimits::VecLimits(index_t size) :     //Construct from size, use default lower bound.
  Low(LOW),
  High(LOW+size-1)
  {
    CHECK;
  }

inline VecLimits::VecLimits(index_t low, index_t high) : //COnstruct from lower and upper bounds.
  Low(low),
  High(high)
  {
    CHECK;
  }

inline VecLimits::~VecLimits() {}
#undef  CHECK

//-----------------------------------------------------------------------------
//
//  Static size calculator
//
inline index_t VecLimits::size(index_t low, index_t high)
{
  return (index_t)(high-low+1);
}

inline index_t VecLimits::size() const
{
  return size(Low,High);
}

//-----------------------------------------------------------------------------
//
//  Allow == operators for comparison (what else?).
//
inline bool VecLimits::operator==(const VecLimits& lim) const
{
  return (Low==lim.Low)&&(High==lim.High);
}

inline bool VecLimits::operator!=(const VecLimits& l) const
{
  return !((*this)==l);
}

inline index_t VecLimits::Offset(index_t i) const
{
  return i-Low;
}

inline std::ostream& operator<<(std::ostream& os,const VecLimits& v)
{
  return v.Write(os);
}

inline std::istream& operator>>(std::istream& is,VecLimits& v)
{
  return v.Read (is);
}

#endif //_veclimit_H_
