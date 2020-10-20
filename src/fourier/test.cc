#include <math.h>
#include <istd::ostream.h>
#include <assert.h>
#include "fft.h"

int main()
{
  int n=2*2*2*2*2*2*2*2*2*2*2*2*2;
  int ni=(int)floor((log(n)/log(2)+1)/2) +2;
  float* wa  =new float[2*n+15];
  int*   ifac=new int  [ni];
  float* r   =new float[n];
  float* rsave=new float[n];
  for (int i=0;i<=n+14;i++) wa[i]=-99;
  for (int i=0;i<ni;i++) ifac[i]=-99;
  for (int i=0;i<n;i++) 
    {
      double i1= (i>=n/2) ? n-i : i;
      double sigma=sqrt(n/10);
      r[i]=sigma/sqrt(M_PI)*exp(-(i1/sigma)*(i1/sigma));
      rsave[i]=r[i];
    }
  drfti1(n,wa+n,ifac);
  drftf1(n,r,wa,wa+n,ifac);
//	for (int i=0;i<=n;i++) std::cout << i << " " << r[i] << std::endl;
	
  float* rt   =new float[n];
  float* ct   =new float[n];
  for (int i=0;i<n;i++) rt[i]=-99;
  for (int i=0;i<n;i++) ct[i]=-99;
  rt[0]=r[0];
  ct[0]=0.0;
  for (int i=1;i<n/2;i++)
  {
    rt[i]=r[2*i-1];
    ct[i]=r[2*i];
    rt[n-i]=rt[i];
    ct[n-i]=ct[i];
  }
  rt[n/2]=r[n-1];
  ct[n/2]=0.0;
  for (int i=0;i<n;i++) std::cout << i << " " << rsave[i] << " " << rt[i] << " " << ct[i] << std::endl;
  std::cout << std::endl;
  std::cout << std::endl;
  return 0;
}
