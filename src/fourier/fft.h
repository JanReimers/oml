// File: fft.h  service routines for the fft package.

void drfti1(int n, float *wa, int *ifac);
void drftf1(int n,float *c,float *ch,float *wa,int *ifac);
void dradf2(int ido,int l1,float *cc,float *ch,float *wa1);
void dradf4(int ido,int l1,float *cc,float *ch,float *wa1,float *wa2,float *wa3);
void dradfg(int ido,int ip,int l1,int idl1,float *cc,float *c1,float *c2,float *ch,float *ch2,float *wa);
