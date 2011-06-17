// File: UT-Random-double.cc  Test distribution of random numbers.

// Copyright (1994-2003), Jan N. Reimers

#include "oml/UnitTests/UnitTest.H"
#include "oml/ran250.h"
#include "oml/histogram.h"
#include "oml/stopw.h"
#include <fstream>
#include <iomanip>

int TestRandomDouble()
{
  bool pass=true;
  const char* Class="Random number generator";
  StartClass(Class);
  StreamableObject::SetToPretty();

  int Nbin=64*64,Nhits=1000;
  Histogram<double> hint(Nbin,0,1.0/(Nbin-1));
  StopWatch sw;
  sw.Start();
  for (int i=0;i<Nbin;i++)
  {
    for (int k=0;k<Nhits;k++) hint.UpDate(OMLRandPos<double>());
  }
  sw.Stop();
  std::cout << sw.GetTime()/(double)(Nbin*Nhits)*1000000.0 << "(ms) per random number." << std::endl;



  Array<double> sigma(hint.GetBins().size());
  for (int i=0;i<sigma.size();i++) sigma[i]=hint.GetBins()[i];
  sigma-=(double)Nhits;
  sigma/=sqrt(Nhits);
  const Array<double>& ss(sigma);
  {
    std::ofstream rawout("raw.dat");
    const Array<int>& sm(hint.GetBins());
    for (int i=0;i<Nbin;i++)
    {
      rawout << i << " " << sm[i] << " " << ss[i] << std::endl;
    }
  }

 double min=-6.0,max=6.0;
 double dsig=1.0/sqrt(Nhits)*4;
  Histogram<double> hsigma(static_cast<int>((max-min)/dsig)+1,min,dsig);
  Array<double>::const_iterator b=sigma.begin();
  for (int i=0;b!=sigma.end();b++,i++)
  {
    if (*b>=min && *b <=max)
      hsigma.UpDate(*b);
    else
      if (i>0 && i<Nbin-1)
      {
        std::cerr << "Large deviation, sigma=" << *b << " bin=" << i << std::endl;
	pass=false;
      }
  }

  std::cout << " x      h(x)    exact(x) " << std::endl;
  const Array<int>& sm(hsigma.GetBins());
  for (int i=0;i<hsigma.GetNumBins();i++)
  {
    double x=min+i*dsig;
    std::cout << x << " " << sm[i] << " " << dsig*Nbin*log(2.0)/sqrt(M_PI)*exp(-x*x*log(2.0)*log(2.0)) << std::endl;
  }


  return pass ? 0 : -1;
}
