#include <cmath>
#include <iostream>
#include <cassert>
#include "fft.h"

// This code comes from netlib.

void fourn(float* data,int* nn,int ndim,int isign);

int main()
{
	int n=2*2*2*2*2;
	float* r   =new float[2*n+1];
	int*   nn  =new int[2];
	nn[1]=n;
	r[0]=-99;
	for (int i=1;i<2*n+1;i+=2) 
	{
		float del=(i-1)/2-n/2.0;
		r[i]=1.0/(2*sqrt(M_PI))*exp(-del*del/4.0);
		r[i+1]=0;
	}
	for (int i=1;i<2*n+1;i+=2) std::cout << (i-1)/2 << " " << r[i] << " " << r[i+1] << std::endl;
  fourn(r,nn,1,1);
	for (int i=1;i<2*n+1;i+=2) std::cout << (i-1)/2 << " " << r[i] << " " << r[i+1] << std::endl;
	std::cout << std::endl;
	std::cout << std::endl;
	return 0;
}
