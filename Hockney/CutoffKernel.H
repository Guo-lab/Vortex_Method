#ifndef _CUTOFFKERNEL_H_
#define _CUTOFFKERNEL_H_

#include "ConvKernel.H"

/**
 * @brief Cutoff kernel implementation for convolutions
 *  Implements a specific type of convolution kernel with a smooth cutoff based on a delta parameter.
 */
class CutoffKernel : public ConvKernel {
  public:
    // Default constructor.
    CutoffKernel() = default;

    /// @brief Destructor
    virtual ~CutoffKernel() = default;

    /// @brief Constructor with parameters
    CutoffKernel(double a_h, double a_delta) : m_h(a_h), m_delta(a_delta) {};

    /**
     * @brief Generate the kernel for convolution
     *
     * @param a_kerArray Output kernel array (pre-allocated)
     * @param a_h Grid spacing
     */
    virtual void getKernel(RectMDArray<std::complex<double>> &a_kerArray, const double &a_h) override;

  private:
    double m_h = 0.0;     ///< Grid spacing
    double m_delta = 0.0; ///< Cutoff radius
};

#endif // _CUTOFFKERNEL_H_
