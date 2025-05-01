#ifndef _FFT1DRecursive_H_
#define _FFT1DRecursive_H_
#include <cmath>
#include <complex>
#include <vector>
#include <cstdio>
#include <iostream>
#include "PowerItoI.H"
#include "FFT1D.H"
class FFT1DRecursive:public FFT1D
{
public:
  FFT1DRecursive(const unsigned int& a_M);
  virtual ~FFT1DRecursive();
  virtual void forwardFFTCC(vector<complex<double> > & a_fHat,
			const vector<complex<double> > & a_f) const;
  virtual void inverseFFTCC(vector<complex<double> > & a_f, 
			const vector<complex<double> > & a_fHat) const;
};
#endif 
