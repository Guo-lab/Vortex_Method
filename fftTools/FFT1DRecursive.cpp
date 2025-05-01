#include <cmath>
#include <complex>
#include <vector>
#include <iostream>
#include <cassert>
#include <stdlib.h>
#include "PowerItoI.H"
#include "FFT1D.H"
#include "FFT1DRecursive.H"

FFT1DRecursive::FFT1DRecursive(const unsigned int& a_M):FFT1D(a_M){};
FFT1DRecursive::~FFT1DRecursive(){};

// Here's a helper function, to avoid "copy&paste" coding. Bool sets forward/backward, just as in FFT1DBRI
void FFTRecursiveHelper(vector<complex<double> >& a_fHat, const vector<complex<double> > & a_f , bool a_isForward, int a_M)
{
  int N = Power(2,a_M);
  assert(N == a_f.size());
  
  int halfSize = N/2;
  assert(halfSize * 2 == a_f.size() );
  vector<complex<double> > fEven,fOdd,fHatEven,fHatOdd;
  fEven.assign(halfSize,0);
  fOdd.assign(halfSize,0);
  fHatEven.assign(halfSize,0);
  fHatOdd.assign(halfSize,0);
  
  for(int i = 0 ; i< N/2 ; i++)
  { 
    fEven[i] = a_f[2*i];
    fOdd[i] = a_f[2*i+1]; 
  }
  if(N > 2)
  {
    FFTRecursiveHelper(fHatEven,fEven,a_isForward,a_M-1);
    FFTRecursiveHelper(fHatOdd,fOdd,a_isForward,a_M-1);
  }
  else if(N == 2)
  {
    assert( fHatEven.size() == 1);
    fHatEven[0] = fEven[0];
    fHatOdd[0] = fOdd[0];
  }
  else{ std::cout << "What the f*ck man. How is N != 2 and < 2?" << std::endl; abort();}


  double length = a_f.size();
  int sign;
  if (a_isForward) {sign =-1;} else {sign=1;}

  complex<double> z0(cos(2*M_PI/length),sign*sin(2*M_PI/length));
  complex<double> zToTheK(1,0);
  for(int k = 0 ; k < N/2 ; k++)
  {
    a_fHat[k] = fHatEven[k] + zToTheK*fHatOdd[k];
    a_fHat[k+N/2] = fHatEven[k] - zToTheK*fHatOdd[k];
    zToTheK *= z0;
  }
}

void FFT1DRecursive::forwardFFTCC(vector<complex<double> >& a_fHat,
			const vector<complex<double> >& a_f) const
{
  assert( a_fHat.size() == m_N);
  bool isForward = true; 
  FFTRecursiveHelper( a_fHat , a_f, isForward, m_M);
}

void FFT1DRecursive::inverseFFTCC(vector<complex<double> >& a_f, const vector<complex<double> > & a_fHat) const
{
  assert(a_f.size() == m_N);
  bool isForward = false;
  FFTRecursiveHelper( a_f, a_fHat, isForward, m_M );
}


