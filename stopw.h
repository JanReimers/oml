// File: stopw.h  StopWatch class from time routines.
#ifndef _stopw_h_
#define _stopw_h_

// Copyright (1994-2003), Jan N. Reimers

#include <time.h>

//#########################################################################
//
//  The StopWatch class is used for timing calculations.   This whole thing
//  is pretty trivial EXCEPT that the number of clock ticks per second
//  seems to be totally UNportable!!!  Every system seems to have a
//  different name.  And on the sun I can't find definition at all!!!!!!
//
class StopWatch
{
  clock_t start,stop;
  double  tics_per_second;
public:
  StopWatch() : start(0), stop(0),tics_per_second
  (
#ifdef __TURBOC__
    CLK_TCK
#endif
#ifdef __GNUC__
    CLOCKS_PER_SEC
#endif
  ) {};
  void   Start  ()       {start=clock();}
  void   Stop   ()       {stop =clock();}
  double GetTime() const {return (double)(stop-start)/tics_per_second;}
};

#endif //_stopw_h_
