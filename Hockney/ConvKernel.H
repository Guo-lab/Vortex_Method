#ifndef _CONVKERNEL_H_
#define _CONVKERNEL_H_

#include <cmath>
#include <complex>
#include <cstdio>
#include <iostream>
#include <vector>

#include "RectMDArray.H"

/**
 * @brief Abstract base class for convolution kernels
 *  Defines the interface for kernels used in the Hockney method for fast convolution in Fourier space.
 */
class ConvKernel {
  public:
    /// Default constructor.
    ConvKernel() {};

    /// Virtual destructor for proper cleanup in derived classes
    virtual ~ConvKernel() = default;

    /**
     * @brief Populate a kernel array (default-constructed RectMDArray) for the convolution operation
     *
     * @param a_srcArray kernel array to be populated
     * @param a_h Grid spacing
     */
    virtual void getKernel(RectMDArray<std::complex<double>> &a_srcArray, const double &a_h) = 0;
};

#endif // _CONVKERNEL_H_
