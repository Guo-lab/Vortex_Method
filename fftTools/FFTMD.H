#ifndef _FFTMD_H_
#define _FFTMD_H_

#include <cassert>
#include <cmath>
#include <iostream>
#include <memory>
#include <vector>

#include "DBox.H"
#include "RectMDArray.H"

#include "FFT1D.H"

/**
 * @brief Multi-dimensional Fast Fourier Transform using 1D FFTs
 *
 * FFTMD implements multi-dimensional FFTs by applying 1D FFTs along each dimension.
 * This class uses a tensor-product approach, where the same 1D FFT is applied
 * along each dimension sequentially.
 */
class FFTMD {
  public:
    /// @brief Default constructor - creates an undefined FFTMD object
    FFTMD() = default;

    /**
     * @brief Constructor with 1D FFT implementation
     *
     * @param a_fft1dPtr Shared pointer to a 1D FFT implementation
     */
    FFTMD(std::shared_ptr<FFT1D> a_fft1dPtr);

    /**
     * @brief Define or redefine an FFTMD object
     *
     * @param a_fft1dPtr Shared pointer to a 1D FFT implementation
     */
    void define(std::shared_ptr<FFT1D> a_fft1dPtr);

    /**
     * @brief Forward FFT on a complex array with standard ordering
     *  Performs forward FFT on input array.
     *
     * @param a_f Complex array to transform (modified in-place)
     */
    void forwardCC(RectMDArray<std::complex<double>> &a_f) const;

    /**
     * @brief Inverse FFT on a complex array with standard ordering
     *
     * @param a_fHat Complex array to transform
     */
    void inverseCC(RectMDArray<std::complex<double>> &a_fHat) const;

    /**
     * @brief Forward FFT on a complex array with centered ordering
     *
     * Performs forward FFT on input array in place, handling centered data
     * where the origin is at the center of the array.
     *
     * @param a_f Complex array to transform (modified in-place)
     */
    void forwardCCcen(RectMDArray<std::complex<double>> &a_f) const;

    /**
     * @brief Inverse FFT on a complex array with centered ordering.
     *
     * @param a_fHat Complex array to transform (modified in-place)
     */
    void inverseCCcen(RectMDArray<std::complex<double>> &a_fHat) const;

    /// Access functions
    inline const int &getN() const { return m_N; };
    inline const int &getM() const { return m_M; };

    /// @brief Check if the FFTMD object is defined
    inline bool isDefined() const { return m_fft1dPtr != nullptr; }

  private:
    int m_N = 0; ///< Number of points in the FFT
    int m_M = 0; ///< Number of dimensions

    std::shared_ptr<FFT1D> m_fft1dPtr; ///< Pointer to the FFT1D object

    /**
     * @brief Allocate vectors for FFT operations
     *
     * @return A pair of vectors for input and output of 1D FFT
     */
    std::pair<std::vector<std::complex<double>>, std::vector<std::complex<double>>> getIOVectors() const {
        assert(isDefined() && "FFTMD::getInputOutputVec FFTMD object not defined");
        return {std::vector<std::complex<double>>(m_N), std::vector<std::complex<double>>(m_N)};
    }

    /**
     * @brief Helper function to create a base DBox for FFT operations
     *  Iterates along a dimension and sets the low and high bounds.
     *
     * @param a_dir Dimension index
     * @param a_N Number of points in the FFT
     * @param a_centered Boolean flag for centered data
     */
    DBox createBaseDBox(int a_dir, bool a_centered) const;

    // Generic data gathering function that works for both standard and centered patterns
    void gatherData(std::vector<std::complex<double>> &a_f1d,
                    const RectMDArray<std::complex<double>> &a_array, const Point &a_pt, const Point &a_edir,
                    bool a_centered) const;

    // Generic data assignment function that works for both standard and centered patterns
    void assignData(RectMDArray<std::complex<double>> &a_array,
                    const std::vector<std::complex<double>> &a_result, const Point &a_pt, const Point &a_edir,
                    bool a_centered) const;
};

#endif
