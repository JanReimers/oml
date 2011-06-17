// File: Angle.H  Overload some comparisons for angles.

#include "Angle.H"
#include <iostream>

const double RingAngle::Pi=M_PI;

std::ostream& operator<<(std::ostream& os,const RingAngle& a)
{
  return os << a.itsValue*180.0/M_PI << "deg";
}
