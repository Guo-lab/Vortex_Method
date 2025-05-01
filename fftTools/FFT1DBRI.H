#ifndef _FFT1DBRI_H_
#define _FFT1DBRI_H_
#include <cmath>
#include <complex>
#include <vector>
#include <cstdio>
#include <iostream>
#include "PowerItoI.H"
#include "FFT1D.H"
class FFT1DBRI:public FFT1D
{
public:  
  FFT1DBRI(const unsigned int& a_M);
  virtual ~FFT1DBRI();
  virtual void forwardFFTCC(vector<complex<double> > & a_fHat, 
                          const vector<complex<double> >& f) const;
  virtual void inverseFFTCC(vector<complex<double> > & a_f, 
                          const vector<complex<double> > & a_fHat) const;
};
#endif
