// File: Angle.H  Overload some comparisons for angles.
#ifndef _Angle_H_
#define _Angle_H_

#include <cmath>
#include <iosfwd>

//
// Imlement ring algebra for angles.
//
class RingAngle
{
public:
  RingAngle() : itsValue(0.0) {};
  RingAngle(double a) : itsValue(a) {};
  const RingAngle& operator=(double a) {itsValue=a;;return *this;}
   
  friend RingAngle operator+ (const RingAngle& a, const RingAngle& b) {return RingAngle(a.itsValue+b.itsValue);}
  friend RingAngle operator- (const RingAngle& a, const RingAngle& b) {return RingAngle(a.itsValue-b.itsValue);}
  friend RingAngle operator* (const RingAngle& a, double b      ) {return RingAngle(a.itsValue*b);}
  friend RingAngle operator/ (const RingAngle& a, double b      ) {return RingAngle(a.itsValue/b);}
  friend bool  operator==(const RingAngle& a, const RingAngle& b) {return fabs((a-b).Normalize())<1e-8;}
  friend bool  operator!=(const RingAngle& a, const RingAngle& b) {return !(a==b);}
  friend bool  operator>=(const RingAngle& a, const RingAngle& b) {return (a-b).Normalize() >=0.0;}
  friend bool  operator<=(const RingAngle& a, const RingAngle& b) {return (b-a).Normalize() >=0.0;}
  friend bool  operator< (const RingAngle& a, const RingAngle& b) {return !(a>=b);}
  friend bool  operator> (const RingAngle& a, const RingAngle& b) {return !(a<=b);}
   
  friend double sin(const RingAngle& a) {return sin(a.itsValue);}
  friend double cos(const RingAngle& a) {return cos(a.itsValue);}
  friend double tan(const RingAngle& a) {return tan(a.itsValue);}
   
  friend std::ostream& operator<<(std::ostream& os,const RingAngle& a);

  double Get() const {return itsValue;}
  static const double Pi;

private:
  double Normalize() //Put in range 0 to 2*Pi.
  {
    itsValue=fmod(itsValue,2*Pi);
    if (itsValue<=-Pi) itsValue+=2*Pi;
    if (itsValue>  Pi) itsValue-=2*Pi;
    return itsValue;
  }
  double itsValue;
};
#endif //_Angle_H_
