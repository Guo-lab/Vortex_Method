#include <chrono>
#include <cmath>
#include <complex>
#include <cstdio>
#include <cstdlib>
#include <iostream>
#include <memory>
#include <string>
#include <vector>

using namespace std;

// #include "FFT1DBRI.H"
// #include "FFT1DRecursive.H"
#include "FFT1D.H"
#include "FFTW1D.H"
#include "PowerItoI.H"

double read_timer() {
    static auto start = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double> elapsed = std::chrono::high_resolution_clock::now() - start;
    return elapsed.count();
}

double test1(shared_ptr<FFT1D> a_fftPtr, double &a_seconds) {
    int N = a_fftPtr->getN();
    vector<complex<double>> f(N);
    vector<double> fSave(N);
    vector<complex<double>> fHat(N);
    double h = 1. / N;

    for (int j = 1; j < N / 2; j++) {
        double x = (j * h - .5) * 32;
        f[j] = complex<double>(exp(-pow(x, 2)), 0.);
        fSave[j] = real(f[j]);
    }

    double time = read_timer();
    a_fftPtr->forwardFFTCC(fHat, f);
    a_fftPtr->inverseFFTCC(f, fHat);
    a_seconds = a_seconds + read_timer() - time;

    double maxerr = 0.;

    for (int j = 0; j < N; j++) {
        if (fabs(real(f[j]) / N - fSave[j]) > maxerr) {
            maxerr = fabs(real(f[j]) / N - fSave[j]);
        }
    }

    return maxerr;
}
int test2(shared_ptr<FFT1D> a_fftPtr, double &a_seconds, int a_mode) {
    int N = a_fftPtr->getN();
    vector<complex<double>> f(N);
    vector<complex<double>> fHat(N);
    double h = 1. / N;
    complex<double> z(cos(2 * M_PI * h * a_mode), sin(2 * M_PI * h * a_mode));
    complex<double> zToTheJ(1., 0.);

    for (int j = 0; j < N; j++) {
        f[j] = zToTheJ;
        zToTheJ *= z;
    }
    double time = read_timer();
    a_fftPtr->forwardFFTCC(fHat, f);
    a_seconds = a_seconds + read_timer() - time;
    int numModes = 0;
    int maxMode = -1;

    for (int j = 0; j < N; j++) {
        if (fabs(real(fHat[j])) / N > 1.e-08) {
            numModes += 1;
            maxMode = j;
        }
    }

    if (numModes > 1) {
        return -1;
    } else {
        return maxMode;
    }
}
int main(int argc, char *argv[]) {
    int M;
    int inputMode;
    string fft_string;

    cout << "input log_2(N), N = number of points" << endl;
    cin >> M;
    cout << "input test mode: recommend something small, but in any case < N" << endl;
    cin >> inputMode;
    cout << "input FFT being tested: Recursive, BRI, FFTW" << endl;
    cin >> fft_string;
    shared_ptr<FFT1D> p_fft;

    // if (fft_string == "Recursive") {
    //     shared_ptr<FFT1DRecursive> p_fft1dR = shared_ptr<FFT1DRecursive>(new FFT1DRecursive(M));
    //     p_fft = dynamic_pointer_cast<FFT1D>(p_fft1dR);
    // } else if (fft_string == "BRI") {
    //     shared_ptr<FFT1DBRI> p_fft1dBRI = shared_ptr<FFT1DBRI>(new FFT1DBRI(M));
    //     p_fft = dynamic_pointer_cast<FFT1D>(p_fft1dBRI);
    // }

    if (fft_string == "FFTW") {
        p_fft = dynamic_pointer_cast<FFT1D>(shared_ptr<FFTW1D>(new FFTW1D(M)));
    } else {
        cout
            << "invalid input - should use BRI, Recursive or FFTW as name for FFT implementation to be tested"
            << endl;
        abort();
    }
    double time = 0.;

    double error = test1(p_fft, time);
    int mode = test2(p_fft, time, inputMode);
    cout << fft_string << ": test 1: error in Gaussian  = " << error << endl;
    cout << fft_string << ": test 2: reproducing input Fourier mode " << inputMode << " , output mode "
         << mode << endl;
    cout << "The input mode number and the output mode number should match" << endl;
    cout
        << "if multiple modes have non-roundoff amplitude, this is an error, and the output mode is set to -1"
        << endl;
    cout << "time spent in calls to fft = " << time << endl;
};
