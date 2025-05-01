#include <cmath>
#include <complex>
#include <cstdio>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <vector>

#include "DBox.H"
#include "FFT1D.H"
#include "FFTMD.H"
#include "FFTW1D.H"

#include "PowerItoI.H"
#include "RectMDArray.H"
#include "WriteRectMDArray.H"

#include "ConvKernel.H"
#include "CutoffKernel.H"
#include "Hockney.H"





/**
 * @brief Test function to check the initialization of Hockney
 */
void testHockneyInitialization() {
    double h = 0.1;
    double delta = 0.2;
    auto cutoffKernel = std::make_shared<CutoffKernel>(h, delta);
    std::shared_ptr<ConvKernel> convkerptr = dynamic_pointer_cast<ConvKernel>(cutoffKernel);

    Hockney hockney(convkerptr, h, 20); // Grid size 2^20 = 1024 * 1024

    assert(hockney.isDefined());
    assert(hockney.getM() == 20);
    assert(hockney.getN() == 1024 * 1024);
    assert(hockney.getH() == h);

    std::cout << std::endl;
    std::cout << "Testing Hockney initialization..." << std::endl;
    std::cout << "Test initialization passed!" << std::endl;
    std::cout << "√" << std::endl;
}

/**
 * @brief Test Hockney solver with a single point source to verify Green's function behavior.
 * Verifies that G * rho matches the theoretical Green's function: ln(r)/(4pi) for r > delta.
 */
void testSinglePointSourceGreenFunction() {
    std::cout << std::endl;
    std::cout << "Testing Hockney with a single point source Green's function..." << std::endl;

    constexpr int M = 7;                       // 2D case
    constexpr int N = 1 << M;                  // Grid size: 64x64
    constexpr double h = 1.0 / N;              // Grid spacing
    constexpr double delta = 0.125;            // Kernel cutoff radius
    constexpr double maxRelativeError = 0.015; // 1.5% error tolerance

    static_assert(DIM == 2, "This test is designed for 2D.");

    // Define the domain
    DBox domain(getZeros(), getOnes() * N);
    RectMDArray<double> rho(domain);
    rho.setVal(0.0);

    // Place a unit charge at the center
    Point center;
    for (int dir = 0; dir < DIM; dir++) center[dir] = N / 2;
    rho[center] = 10000.0;

    // Set up the Hockney solver
    auto cutoffKernel = std::make_shared<CutoffKernel>(h, delta);
    std::shared_ptr<ConvKernel> convkerptr = dynamic_pointer_cast<ConvKernel>(cutoffKernel);
    Hockney hockney(convkerptr, h, M);

    // Perform convolution
    hockney.convolve(rho);

    // Define test points (relative to center) to check the potential
    struct TestPoint {
        int offset_x;
        int offset_y;
    };

    // =========== Test the Green's function with Proportional ratios ================
    std::vector<std::pair<int, int>> offsets = {{16, 0}, {32, 0}, {48, 0}};
    std::vector<Point> testPoints_1;
    for (auto [dx, dy] : offsets) {
        Point p = center;
        p[0] += dx;
        p[1] += dy;
        if (domain.contains(p)) testPoints_1.push_back(p);
    }
    if (testPoints_1.size() == 3) {
        double phi_C = rho[testPoints_1[0]];
        double phi_A = rho[testPoints_1[1]];
        double phi_B = rho[testPoints_1[2]];

        double delta_phi_AC = phi_A - phi_C;
        double delta_phi_BC = phi_B - phi_C;

        double numerical_ratio = delta_phi_AC / delta_phi_BC;
        double theoretical_ratio = 0.63093; // ln(32/16) / ln(48/16) = ln 2 / ln 3

#ifdef DEBUG
        std::cout << "Numerical ratio: " << numerical_ratio << std::endl;
        std::cout << "Theoretical ratio: " << theoretical_ratio << std::endl;
#endif

        if (std::abs(numerical_ratio - theoretical_ratio) < 1e-4)
            std::cout << "Solver behavior is correct! Points: " << testPoints_1[0] << ", " << testPoints_1[1]
                      << " and " << testPoints_1[2] << std::endl;
        else
            assert(false && "Solver behavior is incorrect! There may be an issue");
    }



    // =========== Test the Green's function with equality ================
    auto theoretical_actual_value = [h](double r) { return std::log(r) / (2.0 * M_PI); };
    
    std::vector<std::pair<TestPoint, TestPoint>> testPoints_2 = {
        {{28, 0}, {-28, 0}},    // x-axis symmetry
        {{0, 28}, {0, -28}},    // y-axis symmetry
        {{26, 26}, {-26, -26}}, // diagonal symmetry
        {{28, 0}, {0, 28}}      // cross symmetry
    };
    auto distance = [h](int dx, int dy) { return std::sqrt(dx * dx + dy * dy) * h; };
    for (const auto &pair : testPoints_2) {
        Point p1 = center, p2 = center;
        auto [offset1, offset2] = pair;
        p1[0] += offset1.offset_x;
        p1[1] += offset1.offset_y;
        p2[0] += offset2.offset_x;
        p2[1] += offset2.offset_y;

        if (!domain.contains(p1) || !domain.contains(p2)) continue;

        double r1 = distance(offset1.offset_x, offset1.offset_y); // compute distance from center
        double r2 = distance(offset2.offset_x, offset2.offset_y); // compute distance from center
        if (r1 <= delta && r2 <= delta) {
            std::cout << "Center of the kernel is at (" << center[0] << ", " << center[1] << ")" << std::endl;
            std::cout << "Skipping point (" << p1[0] << ", " << p1[1] << ") and point (" << p2[0] << ", "
                      << p2[1] << ") due to cutoff radius " << delta << std::endl;
            std::cout << "However, the current distance is " << r1 << " and " << r2 << std::endl;
            continue;
        }

        double phi1 = rho[p1];
        double phi2 = rho[p2];
        double difference = std::abs(phi1 - phi2);
        double tolerance = 1e-6;
        if (difference < tolerance)
            std::cout << "Potentials at (" << p1[0] << ", " << p1[1] << "); (" << p2[0] << ", " << p2[1]
                      << ") are equal within tolerance.\n";
        else
            assert(false && "Potentials (Symmetric) differ!\n");

        MDWrite("rho_after_convolve_single_src", rho);
    }

    std::cout << "Single point source Green's function test passed!" << std::endl;
    std::cout << "√" << std::endl;
}



void hockneyTest(const int &a_M) {
    int N = Power(2, a_M);
    double h = 1. / (N);
    double r0 = .25;

    int low[DIM], high[DIM];
    for (int dir = 0; dir < DIM; dir++) {
        low[dir] = 0;
        high[dir] = N;
    }
    DBox b(low, high);
    for (int dir = 0; dir < DIM; dir++) {
        low[dir] = 0;
        high[dir] = 2 * N;
    }
    DBox bDouble(low, high);

    // the source term
    RectMDArray<double> rho(b), rhoSave(b);
    rho.setVal(0.);
    rhoSave.setVal(0.);

    Point center;
    for (int dir = 0; dir < DIM; dir++) center[dir] = N / 2;

    rho[center] = 1.0 / std::pow(h, DIM); // normalized unit charge at center of the grid
    for (Point pt = low; b.notDone(pt); b.increment(pt)) {
        double r = std::sqrt(std::pow(pt[0] * h - .5, 2) + std::pow(pt[1] * h - .5, 2)) / r0;
        if (r < 1.) {
            //    .  .  .  .  .  .  .  .
            //    .  .  .  *  *  .  .  .
            //    .  .  *  *  *  *  .  .
            //    .  *  *  *  *  *  *  .
            //    .  .  *  *  *  *  .  .
            //    .  .  .  *  *  .  .  .
            //    .  .  .  .  .  .  .  .
            //    .  .  .  .  .  .  .  .
            rho[pt] = std::pow(std::cos(M_PI * r / 2), 6);
            rhoSave[pt] = rho[pt];
        }
    }

    double delta = std::pow(h, .75);
    std::cout << "cutoff radius in cells = " << delta / h << std::endl;

    std::shared_ptr<CutoffKernel> cutkptr = std::make_shared<CutoffKernel>(h, delta);
    std::shared_ptr<ConvKernel> convkerptr = dynamic_pointer_cast<ConvKernel>(cutkptr);
    Hockney hockney(convkerptr, h, a_M);

    MDWrite("rho", rho);
    hockney.convolve(rho);
    MDWrite("convRho", rho);
    std::cout << "max value of G * rho " << rho[center] << std::endl;

    DBox bl = b.grow(-2);
    RectMDArray<double> LRho(bl);
    Point unit0 = getUnitv(0);
    Point unit1 = getUnitv(1);

    for (Point pt = bl.getLowCorner(); bl.notDone(pt); bl.increment(pt)) {
        double Laplacian_ninePointStencil =
            (-20 * rho[pt] + (4 * (rho[pt + unit0] + rho[pt - unit0] + rho[pt + unit1] + rho[pt - unit1])) +
             rho[pt + unit0 + unit1] + rho[pt - unit0 + unit1] + rho[pt + unit0 - unit1] +
             rho[pt - unit0 - unit1]) /
            (6 * h * h);
        double averageFiltering =
            rhoSave[pt] + (-4 * rhoSave[pt] + rhoSave[pt + unit0] + rhoSave[pt - unit0] +
                           rhoSave[pt + unit1] + rhoSave[pt - unit1]) /
                              12;
        LRho[pt] = Laplacian_ninePointStencil - averageFiltering;
    }

    MDWrite("residual", LRho);
}


void testProvidedNormTest(int M) {
    std::cout << std::endl;
    std::cout << "Testing Hockney with M since 4 to " << M << std::endl;
    for (int i = 4; i <= M; i++) {

        hockneyTest(i);
#ifdef DEBUG
        std::cout << "Hockney test with M = " << i << " passed." << std::endl;
#endif
    }
    std::cout << "√" << std::endl;
}





/**
 * @brief Test Hockney solver with two point sources to verify
 * wheter the potential at the middle of the two sources (+/-) is zero.
 */
void testTwoPointSources() {
    std::cout << std::endl << "Testing Hockney with two point sources..." << std::endl;

    constexpr int M = 7; // 2D grid, 128x128
    constexpr int N = 1 << M;
    constexpr double h = 1.0 / N;
    constexpr double delta = 0.125;
    constexpr double maxRelativeError = 0.02; // 2% error tolerance

    static_assert(DIM == 2, "This test is designed for 2D.");

    // Define the domain
    DBox domain(getZeros(), getOnes() * N);
    RectMDArray<double> rho(domain);
    rho.setVal(0.0);

    // Place two unit charges at different points
    Point source1, source2;
    for (int dir = 0; dir < DIM; dir++) source1[dir] = N / 4;     // e.g., (32, 32) for N=128
    for (int dir = 0; dir < DIM; dir++) source2[dir] = 3 * N / 4; // e.g., (96, 96) for N=128
    rho[source1] = 10000.0;
    rho[source2] = -10000.0;


    // Set up the Hockney solver
    auto cutoffKernel = std::make_shared<CutoffKernel>(h, delta);
    std::shared_ptr<ConvKernel> convkerptr = dynamic_pointer_cast<ConvKernel>(cutoffKernel);
    Hockney hockney(convkerptr, h, M);

    // Perform convolution
    hockney.convolve(rho);
    Point cancelout;
    for (int dir = 0; dir < DIM; dir++) cancelout[dir] = N / 2; // e.g., (64, 64) for N=128

    double tolerance = 1e-6;
    if (rho[cancelout] < tolerance) {
        std::cout << "The positive (" << source1[0] << ", " << source1[1] << ") and negative (" << source2[0]
                  << ", " << source2[1] << ") charges cancel each other out at the center (" << cancelout[0]
                  << ", " << cancelout[1] << ").\n";
        std::cout << "Two point sources test passed!" << std::endl;
    } else {
        std::cout << "Two point sources test failed." << std::endl;
        assert(false && "Two point sources test failed!");
    }

    MDWrite("rho_after_convolve_pos_neg_src", rho);

    std::cout << "√" << std::endl;
}





int main(int argc, char *argv[]) {
    std::cout << "argc = " << argc << std::endl;
    if (argc < 2) {
        std::cout << "Usage: " << argv[0] << " <M>" << std::endl;
        return 1;
    }

    int M;
    sscanf(argv[1], "%d", &M);

    testHockneyInitialization();

    testSinglePointSourceGreenFunction();
    testTwoPointSources();

    testProvidedNormTest(M);

    return 0;
}
