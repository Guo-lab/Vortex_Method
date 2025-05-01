#include "WriteRectMDArray.H"

#include "Hockney.H"
#include "PowerItoI.H"





Hockney::Hockney(std::shared_ptr<ConvKernel> &a_kerPtr, const double &a_h, int a_M) {
    define(a_kerPtr, a_h, a_M);
}

void Hockney::define(std::shared_ptr<ConvKernel> &a_kerPtr, const double &a_h, int a_M) {
    assert(a_kerPtr != nullptr && "Kernel pointer cannot be null");
    assert(a_M > 0.0 && "Grid size power M must be greater than 0");
    assert(a_h > 0.0 && "Grid spacing h must be greater than 0");

    std::shared_ptr<FFTW1D> p_fftw1d = std::make_shared<FFTW1D>(a_M + 1);
    std::shared_ptr<FFT1D> p_fft = std::dynamic_pointer_cast<FFT1D>(p_fftw1d);

    m_fftmd.define(p_fft);
    m_kerPtr = a_kerPtr;
    m_h = a_h;
    m_M = a_M;
    m_N = Power(2, a_M);
    m_isDefined = true;
}





void Hockney::validateInputDomain(const RectMDArray<double> &a_rhs) const {
    assert(m_isDefined && "Hockney solver not properly defined");

    DBox rhsDomain = a_rhs.getDBox();
    Point low = rhsDomain.getLowCorner();
    Point high = rhsDomain.getHighCorner();

    assert(low == getZeros() && "Input domain must start at origin");
    assert(high == getOnes() * m_N && "Input domain dimensions must match initialization N");
}


std::pair<RectMDArray<std::complex<double>>, RectMDArray<std::complex<double>>>
Hockney::prepareComplexArrays(const RectMDArray<double> &a_rhs, const DBox &a_expandedDomain) const {
    RectMDArray<std::complex<double>> rhsComplex(a_expandedDomain);
    RectMDArray<std::complex<double>> kernelComplex(a_expandedDomain);

    std::complex<double> zero(0.0, 0.0);
    rhsComplex.setVal(zero);


    DBox rhsDomain = a_rhs.getDBox();
    for (Point pt = rhsDomain.getLowCorner(); rhsDomain.notDone(pt); rhsDomain.increment(pt)) {
        rhsComplex[pt].real(a_rhs[pt]);
    }


    m_kerPtr->getKernel(kernelComplex, m_h);

    return {rhsComplex, kernelComplex};
}





void Hockney::performFFTOperations(RectMDArray<std::complex<double>> &a_rhsComplex,
                                   RectMDArray<std::complex<double>> &a_kernelComplex) const {
    m_fftmd.forwardCCcen(a_rhsComplex);
    m_fftmd.forwardCCcen(a_kernelComplex);

    DBox domain = a_rhsComplex.getDBox();
    for (Point pt = domain.getLowCorner(); domain.notDone(pt); domain.increment(pt)) {
        a_rhsComplex[pt] *= a_kernelComplex[pt];
    }

    m_fftmd.inverseCCcen(a_rhsComplex);
}





void Hockney::convolve(RectMDArray<double> &a_rhs) {
    validateInputDomain(a_rhs);
    DBox expandedDomain = createExpandedDomain(a_rhs.getDBox());

    // To compute the linear convolution of two signals of size N * N,
    // need to zero-pad them to at least (2N − 1) * (2N − 1).
    auto [rhsDouble, kernel] = prepareComplexArrays(a_rhs, expandedDomain);

    performFFTOperations(rhsDouble, kernel);

    a_rhs.setVal(0.);
    DBox rhsDomain = a_rhs.getDBox();

    // assumes periodic boundary conditions;
    // avoids wrap-around effect
    DBox validRegion(rhsDomain.getLowCorner(), rhsDomain.getHighCorner() - getOnes());

    // A forward FFT followed by an inverse FFT,
    // original values multiplied by N in one dimension.
    // Domain twice the original size in each dimension,
    // so DIM to the power 2.
    double scale = 1. / pow(m_N * 1., DIM * 2) / 4;

    for (Point pt = validRegion.getLowCorner(); validRegion.notDone(pt); validRegion.increment(pt)) {
        a_rhs[pt] = real(rhsDouble[pt]) * scale;
    }
}
