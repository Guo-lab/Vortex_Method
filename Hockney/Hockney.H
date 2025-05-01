#ifndef _HOCKNEY_
#define _HOCKNEY_

#include <array>
#include <assert.h>
#include <cmath>
#include <complex>
#include <cstdio>
#include <iostream>
#include <memory>
#include <stdlib.h>
#include <vector>

#include "ConvKernel.H"

#include "FFT1D.H"
#include "FFTMD.H"
#include "FFTW1D.H"

/**
 * @brief Fast convolution solver using the Hockney method
 *
 * Implements fast convolution via FFTs for solving Poisson-like
 * equations and other convolution problems.
 */
class Hockney {
  public:
    /// Default constructor - creates undefined Hockney object
    Hockney() { m_isDefined = false; };

    /// Constructor - creates Hockney object with specified kernel and grid spacing
    Hockney(std::shared_ptr<ConvKernel> &a_kerPtr, const double &a_h, int a_M);

    /**
     * @brief Define or redefine Hockney object
     *
     * @param a_kerPtr Convolution kernel
     * @param a_h Grid spacing
     * @param a_M Power of 2 for grid size (N = 2^M)
     */
    void define(std::shared_ptr<ConvKernel> &a_kerPtr, const double &a_h, int a_M);

    /**
     * @brief Validate that the input domain meets the requirements
     *
     * @param a_rhs Input array to validate
     * @throws std::invalid_argument If domain does not meet requirements
     */
    void validateInputDomain(const RectMDArray<double> &a_rhs) const;

    /**
     * @brief Create expanded domain for FFT operations
     *
     * @param a_origDomain Original domain
     * @return DBox Expanded domain for FFT operations
     */
    DBox createExpandedDomain(const DBox &a_origDomain) const {
        return DBox(a_origDomain.getHighCorner() * (-1), a_origDomain.getHighCorner());
    }

    /**
     * @brief Prepare complex arrays for FFT operations
     *
     * @param a_rhs Original real array
     * @param a_expandedDomain Domain for complex arrays
     * @return std::pair<RectMDArray<std::complex<double>>, RectMDArray<std::complex<double>>>
     *         Pair of complex arrays (input and kernel)
     */
    std::pair<RectMDArray<std::complex<double>>, RectMDArray<std::complex<double>>>
    prepareComplexArrays(const RectMDArray<double> &a_rhs, const DBox &a_expandedDomain) const;

    /**
     * @brief Perform FFT operations and frequency domain multiplication
     *
     * @param a_rhsComplex Input/output complex array
     * @param a_kernelComplex Kernel complex array
     */
    void performFFTOperations(RectMDArray<std::complex<double>> &a_rhsComplex,
                              RectMDArray<std::complex<double>> &a_kernelComplex) const;

    /**
     * @brief Perform convolution using the Hockney method
     *
     * i.      Extend the domain to handle periodic boundary conditions
     *  ii.     Convert to complex values for FFT operations
     *  iii.    Transform to frequency domain
     *  iv.     Apply convolution (multiplication in frequency domain)
     *  v.      Transform back to spatial domain
     *  vi.     Extract and scale the result
     *
     * @param a_rhs Right-hand side Input array to be convolved (modified in-place)
     */
    void convolve(RectMDArray<double> &a_rhs);

    /// @brief Check if Hockney object is properly defined
    bool isDefined() const { return m_isDefined; }

    int getM() const { return m_M; }    ///< Get power of 2 for grid size
    int getN() const { return m_N; }    ///< Get grid size (N = 2^M)
    double getH() const { return m_h; } ///< Get grid spacing

    const ConvKernel &getKernel() const { return *m_kerPtr; } ///< Get convolution kernel

    /// Destructor
    ~Hockney() = default;

  protected:
    double m_h = 0.0; ///< Grid spacing
    int m_M;          ///< Power of 2 for grid size
    int m_N;          ///< Grid size (N = 2^M)

    FFTMD m_fftmd;                        ///< FFTMD object for performing FFTs
    std::shared_ptr<ConvKernel> m_kerPtr; ///< Convolution kernel

    bool m_isDefined = false; ///< Flag for proper initialization
};

#endif // _HOCKNEY_
