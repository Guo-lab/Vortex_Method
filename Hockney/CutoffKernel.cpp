#include "CutoffKernel.H"

/**
 * @ref
 * https://math.libretexts.org/Bookshelves/Differential_Equations/Introduction_to_Partial_Differential_Equations_(Herman)/07%3A_Green's_Functions/7.05%3A_Greens_Functions_for_the_2D_Poisson_Equation
 *
 * The cutoff should be applied around the source point, not the origin. In a convolution-based solver,
 * the kernel is defined with a cutoff around its center (e.g., [0,0]),
 * but the convolution shifts this center to each source point’s location,
 * correctly placing the cutoff around the source regardless of its position.
 *
 *
 */
void CutoffKernel::getKernel(RectMDArray<std::complex<double>> &a_kerArray, const double &a_h) {

    DBox domain = a_kerArray.getDBox();
    Point low = domain.getLowCorner();
    Point high = domain.getHighCorner();

    // Initialize kernel array to zero
    std::complex<double> one(0., 0.);
    a_kerArray.setVal(one);

    // Calculate squared delta
    double deltaSquare = std::pow(m_delta, 2);

    // Fill the kernel array based on radial distance
    for (Point pt = low; domain.notDone(pt); domain.increment(pt)) {

        // Calculate squared radius normalized by delta ot the power of 2
        // Computing the squared Euclidean distance of the current point from the origin,
        // and dividing the squared distance.
        double rSq = (std::pow(pt[0] * m_h, 2) + std::pow(pt[1] * m_h, 2)) / deltaSquare;

        if (rSq < 1.0) {

            // Inside cutoff radius - use 7th-like order polynomial approximation for 1/(2pi) lnr
            double r = std::sqrt(rSq);
            double kernelValue =
                (((((2160 * r - 9800) * r + 16464) * r - 11025) * r * r + 2940) * r * r - 739.) / 420.0 /
                (2.0 * M_PI);
            a_kerArray[pt].real(kernelValue);

        } else {

            // Outside cutoff radius - use logarithmic formula in a infinite plane
            a_kerArray[pt].real(std::log(rSq) / (4 * M_PI));
        }
    }
}