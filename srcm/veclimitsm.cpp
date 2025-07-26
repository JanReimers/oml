module;
#include <cstdint>
#include <iostream>
#include <cassert>

export module oml.VecLimits;
import oml.StreamableObject;


export typedef int64_t index_t;


// veclimits.h
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
export class VecLimits
{
 public:
  /*! \name Constructors*/
  //@{
  VecLimits(               ); //!<Limits for a null Vector.
  VecLimits(size_t         ); //!<Construct from size, use default lower bound.
  VecLimits(index_t,index_t); //!<Construct from lower and upper bounds.
  //@}
 ~VecLimits();

  void ReBase(int low);

  static size_t  size(index_t,index_t)      ;
  //! Returns number of elements.
  size_t  size(               ) const;

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

inline VecLimits::VecLimits(size_t size) :     //Construct from size, use default lower bound.
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

inline void VecLimits::ReBase(int low)
{
    High=High-Low+low;
    Low=low;
}

//-----------------------------------------------------------------------------
//
//  Static size calculator
//
inline size_t  VecLimits::size(index_t low, index_t high)
{
  return (index_t)(high-low+1);
}

inline size_t  VecLimits::size() const
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

//-----------------------------------------------------------------------------
//
//  IO
//
std::ostream& VecLimits::Write(std::ostream& os) const
{
  assert(os);
  if(StreamableObject::Binary()) {BinaryWrite(Low,os);BinaryWrite(High,os);}
  if(StreamableObject::Ascii ()) os << Low << " " << High << " ";
  if(StreamableObject::Pretty()) os << "(" << Low << ":" << High << ")";
  assert(os);
  return os;
}

std::istream& VecLimits::Read(std::istream& is)
{
  assert(is);
  if(StreamableObject::Binary()) {BinaryRead(Low,is);BinaryRead(High,is);}
  if(StreamableObject::Ascii ()) is >> Low >> High;
  assert(is);
  return is;
}



export inline std::ostream& operator<<(std::ostream& os,const VecLimits& v)
{
  return v.Write(os);
}

export inline std::istream& operator>>(std::istream& is,VecLimits& v)
{
  return v.Read (is);
}

// Use this to get limits of a tensor product
export VecLimits operator*(const VecLimits& a, const VecLimits& b);


//-----------------------------------------------------------------------------
//
//  Index checking with descriptive erros.
//
bool VecLimits::Check() const
{
  return High+1 >= Low;
}

bool VecLimits::CheckIndex(index_t i) const
{
  bool ok=true;
  if (i< Low )
  {
    std::cout << "Index " << i << " too low in VecLimits, Low=" << Low << std::endl;
    ok=false;
  }
  if (i>High)
  {
    std::cout << "Index " << i << " too high in VecLimits, High=" << High << std::endl;
    ok=false;
  }
#ifdef DEBUG
  assert(ok);
#endif
  return ok;
}

export VecLimits operator*(const VecLimits& a, const VecLimits& b)
{
    assert(a.Low==b.Low);
    int size=a.size()*b.size();
    return VecLimits(a.Low,a.Low+size-1);
}


#undef CHECK
