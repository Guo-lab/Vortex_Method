#include <cassert>
#include <cmath>
#include <cstdio>
#include <iostream>
#include <vector>

#include "CH_Timer.H"
#include "DBox.H"
#include "RectMDArray.H"

#include "FFT1D.H"
#include "FFTMD.H"

FFTMD::FFTMD(std::shared_ptr<FFT1D> a_fft1dPtr) { define(a_fft1dPtr); }

void FFTMD::define(std::shared_ptr<FFT1D> a_fft1dPtr) {
    assert(a_fft1dPtr != nullptr);

    m_fft1dPtr = a_fft1dPtr;
    m_M = m_fft1dPtr->getM();
    m_N = m_fft1dPtr->getN();
}

void FFTMD::forwardCC(RectMDArray<std::complex<double>> &a_f) const {
    assert(isDefined() && "FFTMD::forwardCC FFTMD object not defined");

    auto [f1d, fHat1d] = getIOVectors();

    CH_TIMERS("fftmdForward");
    CH_TIMER("fft1d", t1);

    for (int dir = 0; dir < DIM; dir++) {
        DBox base = createBaseDBox(dir, false);
        Point edir = getUnitv(dir);

        for (Point pt = base.getLowCorner(); base.notDone(pt); base.increment(pt)) {
            gatherData(f1d, a_f, pt, edir, false);
            CH_START(t1);
            m_fft1dPtr->forwardFFTCC(fHat1d, f1d);
            CH_STOP(t1);
            assignData(a_f, fHat1d, pt, edir, false);
        }
    }
}

void FFTMD::inverseCC(RectMDArray<std::complex<double>> &a_fHat) const {
    assert(isDefined() && "FFTMD::inverseCC FFTMD object not defined");

    auto [f1d, fHat1d] = getIOVectors();

    CH_TIMERS("fftmdInverse");
    CH_TIMER("fft1drev", t3);

    for (int dir = 0; dir < DIM; dir++) {
        DBox base = createBaseDBox(dir, false);
        Point edir = getUnitv(dir);

        for (Point pt = base.getLowCorner(); base.notDone(pt); base.increment(pt)) {
            gatherData(fHat1d, a_fHat, pt, edir, false);
            CH_START(t3);
            m_fft1dPtr->inverseFFTCC(f1d, fHat1d);
            CH_STOP(t3);
            assignData(a_fHat, f1d, pt, edir, false);
        }
    }
}

void FFTMD::forwardCCcen(RectMDArray<std::complex<double>> &a_f) const {
    assert(isDefined() && "FFTMD::forwardCCcen FFTMD object not defined");

    auto [f1d, fHat1d] = getIOVectors();

    CH_TIMERS("fftmdForward");
    CH_TIMER("fft1d", t1);

    for (int dir = 0; dir < DIM; dir++) {
        DBox base = createBaseDBox(dir, true);
        Point edir = getUnitv(dir);

        for (Point pt = base.getLowCorner(); base.notDone(pt); base.increment(pt)) {
            gatherData(f1d, a_f, pt, edir, true);
            CH_START(t1);
            m_fft1dPtr->forwardFFTCC(fHat1d, f1d);
            CH_STOP(t1);
            assignData(a_f, fHat1d, pt, edir, true);
        }
    }
}

void FFTMD::inverseCCcen(RectMDArray<std::complex<double>> &a_fHat) const {
    assert(isDefined() && "FFTMD::inverseCCcen FFTMD object not defined");

    auto [f1d, fHat1d] = getIOVectors();

    CH_TIMERS("fftmdInverse");
    CH_TIMER("fft1drev", t3);

    for (int dir = 0; dir < DIM; dir++) {
        DBox base = createBaseDBox(dir, true);
        Point edir = getUnitv(dir);

        for (Point pt = base.getLowCorner(); base.notDone(pt); base.increment(pt)) {
            gatherData(fHat1d, a_fHat, pt, edir, true);

            CH_START(t3);
            m_fft1dPtr->inverseFFTCC(f1d, fHat1d);
            CH_STOP(t3);

            assignData(a_fHat, f1d, pt, edir, true);
        }
    }
}

/// ========================= HELPER FUNCTIONS =========================
DBox FFTMD::createBaseDBox(int a_dir, bool a_centered) const {
    int low[DIM], high[DIM];

    for (int dir2 = 0; dir2 < DIM; dir2++) {
        low[dir2] = a_centered ? -m_N / 2 : 0;
        high[dir2] = a_centered ? m_N / 2 - 1 : m_N - 1;
    }
    high[a_dir] = 0;
    if (a_centered) low[a_dir] = 0;

    Point lowCorner(low), highCorner(high);
    DBox base(lowCorner, highCorner);
    return base;
}

void FFTMD::gatherData(std::vector<std::complex<double>> &a_f1d,
                       const RectMDArray<std::complex<double>> &a_array, const Point &a_pt,
                       const Point &a_edir, bool a_centered) const {
    if (a_centered) {
        for (int l = 0; l < m_N / 2; l++) {
            a_f1d[l] = a_array[a_pt + a_edir * l];
            a_f1d[m_N - l - 1] = a_array[a_pt - a_edir * (l + 1)];
        }
    } else {
        for (int l = 0; l < m_N; l++)
            a_f1d[l] = a_array[a_pt + a_edir * l];
    }
}

void FFTMD::assignData(RectMDArray<std::complex<double>> &a_array,
                       const std::vector<std::complex<double>> &a_result, const Point &a_pt,
                       const Point &a_edir, bool a_centered) const {
    if (a_centered) {
        for (int l = 0; l < m_N / 2; l++) {
            a_array[a_pt + a_edir * l] = a_result[l];
            a_array[a_pt - a_edir * (l + 1)] = a_result[m_N - l - 1];
        }
    } else {
        for (int l = 0; l < m_N; l++)
            a_array[a_pt + a_edir * l] = a_result[l];
    }
}
