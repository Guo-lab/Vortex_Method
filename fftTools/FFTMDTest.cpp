#include <cmath>
#include <complex>
#include <cstdio>
#include <iostream>
#include <memory>
#include <vector>

#include "DBox.H"
#include "RectMDArray.H"

#include "PowerItoI.H"

#include "FFT1D.H"
#include "FFTMD.H"
#include "FFTW1D.H"

// #include "FFT1DBRI.H"

using namespace std;

double test(const FFTMD &a_fftmd) {
    int N = a_fftmd.getN();
    int low[DIM], high[DIM];

    for (int dir = 0; dir < DIM; dir++) {
        low[dir] = 0;
        high[dir] = N - 1;
    }
    DBox b(low, high);

    RectMDArray<complex<double>> f(b);
    RectMDArray<complex<double>> fSave(b);
    RectMDArray<complex<double>> fHat(b);
    double h = 1. / N;

    for (Point pt = b.getLowCorner(); b.notDone(pt); b.increment(pt)) {
        double x = 1.;
        for (int dir = 0; dir < DIM; dir++) {

            double y = (pt[dir] * h - .5);
            x += y * y * 32 * 32;
        }
        f[pt] = complex<double>(exp(-x), 0);
        fSave[pt] = f[pt];
    }
    double maxramp = 0.;
    double maxiamp = 0.;
    double maxamp = 0.;
    Point ptmax = getZeros();

    a_fftmd.forwardCC(f);
    for (Point pt = b.getLowCorner(); b.notDone(pt); b.increment(pt)) {
        if (real(f[pt]) * real(f[pt]) + imag(f[pt]) * imag(f[pt]) > maxamp) {
            maxramp = real(f[pt]);
            maxiamp = imag(f[pt]);
            maxramp = real(f[pt]) * real(f[pt]) + imag(f[pt]) * imag(f[pt]);
            ptmax = pt;
        }
    }

    a_fftmd.inverseCC(f);

    double maxerr = 0.;
    double minerr = 10000.;
    int normalize = Power(N, DIM);

    for (Point pt = b.getLowCorner(); b.notDone(pt); b.increment(pt)) {
        if (fabs(real(f[pt]) / normalize) > maxerr) {
            maxerr = fabs(real(f[pt]) / normalize - real(fSave[pt]));
        }
    }
    return maxerr;
}

int main(int argc, char *argv[]) {
    int M;
    double time;
    string fft_string;
    cout << "input log_2(N), N = number of points" << endl;
    cin >> M;
    cout << "input FFT being tested: BRI or FFTW" << endl;
    cin >> fft_string;

    shared_ptr<FFT1D> p_fft;

    // if (fft_string == "BRI") {
    //     p_fft = dynamic_pointer_cast<FFT1D>(shared_ptr<FFT1DBRI>(new FFT1DBRI(M)));
    // }
    if (fft_string == "FFTW") {
        p_fft = dynamic_pointer_cast<FFT1D>(shared_ptr<FFTW1D>(new FFTW1D(M)));
    }

    if (fft_string != "BRI" && fft_string != "FFTW") {
        cout << "invalid input - should use BRI or FFTW as name for FFT implementation to be tested" << endl;
        abort();
    }

    FFTMD fftmd(p_fft);
    double error = test(fftmd);

    cout << "test: error in Gaussian  = " << error << endl;
};
