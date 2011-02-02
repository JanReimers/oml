// File: histogram.h  inline histogram class.
#ifndef _histogram_h_
#define _histogram_h_

// Copyright (1994-2005), Jan N. Reimers

#include "oml/array.h"

template <class T> class Histogram
{
 public:
  Histogram();  
  Histogram(int nbin, T xmin, T dx);  

  void   UpDate(const T& X);
  bool   InRange(const T& X) const;
  const Array<int>& GetBins() const {return itsBins;}
  Array<T> GetXs() const;
  int GetNumBins() const {return itsNumBins;}

  double GetMoment  (int n)  const;
  double GetAverage() const {return GetMoment(1);}
  double GetSigma  () const 
  {
    double hbar=GetAverage();
    return sqrt((GetMoment(2)-hbar*hbar)*itsNumHits/(itsNumHits-1));
  }
  
 private:

  double GetOverlap (const Histogram&) const;

  VecLimits GetPackedLimits() const;
  VecLimits GetUnpackedLimits() const {return VecLimits(0,itsNumBins-1);}
  void      Pack(); //Compress (ignore many 0 count bins) bin arrays for storage.
  void      UnPack();
   
  T   itsXMin, itsXMax,itsDX,itsScale;
  int itsNumHits, itsNumBins;
  Array<int> itsBins;
};

template <class T> Histogram<T>::Histogram()
  : itsXMin       (0)
  , itsXMax       (0)
  , itsDX         (0)
  , itsScale      (0)
  , itsNumHits    (0)
  , itsNumBins    (0)
  , itsBins       ( ) 
  {}

template <class T> Histogram<T>::Histogram(int nbin, T xmin, T dx) 
  : itsXMin       (xmin)
  , itsXMax       (xmin+nbin*dx)
  , itsDX         (dx)
  , itsScale      (T(1)/itsDX)
  , itsNumHits    (0)
  , itsNumBins    (nbin)
  , itsBins       (nbin)        //C-style subscripting for bins!
{
   assert(itsDX>0);
   assert(nbin>0);
   Fill(itsBins,0);
}

//
//  Only a new histogram with all bins allocated can be updated.
//  Once the histogram has been saved to disk and restored, the
//  bin arrays limits will be compressed. 
//
template <class T> inline void Histogram<T>::UpDate(const T& X)
{
  index_t bin=(index_t)floor((X-itsXMin)*itsScale+0.5);
  if (bin>=0&&bin<itsBins.size())
  {
    itsBins[bin]++;
    itsNumHits++;
  }
}

template <class T> inline bool Histogram<T>::InRange(const T& X) const
{
  return X>=itsXMin && X<itsXMax; //THis is not self consistent with the update algorithm, and therefore useless!
}

template <class T> Array<T> Histogram<T>::GetXs() const
{
  int N=GetNumBins();
  Array<T> ret(N);
  for (int i=0;i<N;i++) ret[i]=itsXMin+i*itsDX;
  return ret;
}

template <class T> double Histogram<T>::GetMoment  (int n)  const
{
  assert(itsBins.size()==itsNumBins);
  T ret(0);
  for (int i=0;i<itsNumBins;i++)
  {
    T X=itsXMin+i*itsDX;
    ret+=pow(X,n)*itsBins[i];
  }
  return ret/itsNumHits;
}


template <class T> double Histogram<T>::GetOverlap(const Histogram& h) const
{
  int low  = Max(0,h.itsBins.size()-1);
  int high = Min(0,h.itsBins.size()-1);
   
  int ret=0;
  for (int i=low;i<=high;i++) ret+=itsBins(i)+h.itsBins(i);
  return 100.0*ret/(itsNumHits+h.itsNumHits); //Percnetage of total hits.
}


template <class T> VecLimits Histogram<T>::GetPackedLimits() const
{
   int low =0;
   int high=itsBins.size()-1;
   for (int i=low ;i<=high;i++) if (itsBins(i)>0) {low =i;break;}
   for (int i=high;i>=low ;i--) if (itsBins(i)>0) {high=i;break;}
   return VecLimits(low,high);
}

/*template <class T> inline void Histogram<T>::Pack()
{
   itsBins.SetLimits(GetPackedLimits(),true);
}

template <class T> void Histogram<T>::UnPack()
{
   VecLimits oldlim=itsBins.GetLimits();
   itsBins.SetLimits(GetUnpackedLimits(),true);

   Array<int   >::Subscriptor sb(itsBins);
   for(subsc_t i=GetUnpackedLimits().Low; i<oldlim.Low; i++) sb[i]=0;
   for(subsc_t i=oldlim.High+1; i<=GetUnpackedLimits().High; i++) sb[i]=0;
   }*/

#endif // _histogram_h_
