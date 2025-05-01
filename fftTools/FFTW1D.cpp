#include "FFTW1D.H"

FFTW1D::FFTW1D(unsigned int a_M) : FFT1D(a_M) {
    m_in.resize(m_N);
    m_out.resize(m_N);
    fftw_complex *in;
    fftw_complex *out;
    in = reinterpret_cast<fftw_complex *>(&(m_in[0]));
    out = reinterpret_cast<fftw_complex *>(&(m_out[0]));

    m_forward = fftw_plan_dft_1d(m_N, in, out, FFTW_FORWARD, FFTW_ESTIMATE);
    m_inverse = fftw_plan_dft_1d(m_N, in, out, FFTW_BACKWARD, FFTW_ESTIMATE);
}

FFTW1D::~FFTW1D() {
    fftw_destroy_plan(m_forward);
    fftw_destroy_plan(m_inverse);
}

void FFTW1D::forwardFFTCC(std::vector<std::complex<double>> &a_fHat,
                          const std::vector<std::complex<double>> &f) const {
    m_in = f;
    fftw_execute(m_forward);
    a_fHat = m_out;
}

void FFTW1D::inverseFFTCC(std::vector<std::complex<double>> &a_f,
                          const std::vector<std::complex<double>> &a_fHat) const {
    m_in = a_fHat;
    fftw_execute(m_inverse);
    a_f = m_out;
}
