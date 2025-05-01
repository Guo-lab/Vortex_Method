#include <cmath>
#include <complex>
#include <vector>
#include <cstdio>
#include <iostream>
#include <cassert>
#include "PowerItoI.H"
#include "FFTCTBRI.H"
#include "FFT1D.H"
#include "FFT1DBRI.H"
FFT1DBRI::FFT1DBRI(const unsigned int& a_M):FFT1D(a_M)
{
};
FFT1DBRI::~FFT1DBRI(){};
void FFT1DBRI::forwardFFTCC(vector<complex<double> > & a_fHat, 
                            const vector<complex<double> >& a_f) const
{
  assert (a_fHat.size()==m_N);
  for (unsigned int j = 0; j < m_N; j++)
    {
      a_fHat[j]= a_f[j];
    }
  bool isForward = true;
  
  FFTCTBRI(a_fHat, m_M, isForward);
}
void FFT1DBRI::inverseFFTCC(vector<complex<double> >& a_f , 
                                  const vector<complex<double> > & a_fHat ) const
{
  assert (a_fHat.size()==m_N);
  for (unsigned int j = 0; j < m_N; j++)
    {
      a_f[j]=a_fHat[j];
    }
  bool isForward = false;
  FFTCTBRI(a_f, m_M, isForward);
}
