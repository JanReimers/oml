#include "oml/vector.h"
#include "oml/imp/iterable.h"
#include "oml/random.h"


#include <iostream>
#include <iomanip>
#include <cassert>

#if DEBUG
#define CHECK Check()
#else
#define CHECK
#endif

//-----------------------------------------------------------------------------
//
//  IO
//
template <class T> std::ostream& Vector<T>::Write(std::ostream& os) const
{
  assert(os);
  if (!this->Pretty())
    os << GetLimits();
   else
  {
    std::streamsize prec=os.precision();
    std::streamsize wid =os.width();
    os << GetLimits();
    os << std::setw(wid) << std::setprecision(prec);
  }
  return ::Write(os,*this);
}

template <class T> std::istream& Vector<T>::Read(std::istream& is)
{
  assert(is);
  VecLimits lim;
  is >> lim;
  if (size()==0)
    SetLimits(lim);
   else
    assert(GetLimits()==lim);

  ::Read(is,*this);
  CHECK;
  assert(is);
  return is;
}

#undef CHECK

