#ifndef _FFTW1D_H_
#define _FFTW1D_H_

#include "fftw3.h"
#include <cmath>
#include <complex>
#include <cstdio>
#include <cstdlib>
#include <iostream>
#include <vector>

#include "FFT1D.H"
#include "PowerItoI.H"

class FFTW1D : public FFT1D {
  public:
    /// @brief Default constructor.
    FFTW1D(unsigned int a_M);

    /// @brief Destructor.
    virtual ~FFTW1D();

    virtual void forwardFFTCC(std::vector<std::complex<double>> &a_fHat,
                              const std::vector<std::complex<double>> &f) const;

    virtual void inverseFFTCC(std::vector<std::complex<double>> &a_f,
                              const std::vector<std::complex<double>> &a_fHat) const;

  protected:
    mutable std::vector<std::complex<double>> m_in; ///< Input data for FFTW
    std::vector<std::complex<double>> m_out;        ///< Output data for FFTW
    fftw_plan m_forward, m_inverse;                 ///< FFTW plans
};

#endif // _FFTW1D_H_
