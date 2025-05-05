#include "CutoffKernel.H"
#include "ParticleSet.H"
#include "ParticleVelocities.H"
#include "RK4.H"
#include "RectMDArray.H"
#include "VisitWriter.H"
#include "WriteRectMDArray.H"
#include <cassert>
#include <cmath>
#include <iostream>
#include <vector>





void outField(ParticleSet &p, int a_coarsenFactor) {
    int coarsenFactor = a_coarsenFactor;

    DBox bx = p.m_box.coarsen(coarsenFactor);
    double h = p.m_dx * coarsenFactor;

    RectMDArray<double> outVort(bx);

    array<int, DIM> ipos;
    array<double, DIM> xpos;
    double weight;

    Point e0 = getUnitv(0);
    Point e1 = getUnitv(1);

    outVort.setVal(0.);
    for (int k = 0; k < p.getSize(); k++) {
        for (int l = 0; l < DIM; l++) {
            double newpos = p.getParticle(k).m_x[l];
            ipos[l] = newpos / h;
            xpos[l] = (newpos - ipos[l] * h) / h;
        }

        Point pt(ipos);

#ifdef DEBUG
        std::cout << "pt = " << pt << std::endl;
#endif
        assert(p.m_box.contains(pt));
#ifdef DEBUG
        std::cout << "pt = " << pt << std::endl;
#endif

        for (int l0 = 0; l0 < DIM; l0++) {
            for (int l1 = 0; l1 < DIM; l1++) {
                outVort[pt + e0 * l0 + e1 * l1] += (1. - xpos[0] + (2 * xpos[0] - 1.) * l0) *
                                                   (1. - xpos[1] + (2 * xpos[1] - 1.) * l1) *
                                                   p.getParticle(k).strength / coarsenFactor / coarsenFactor;
            }
        }
    }

    const char *foo = MDWrite(outVort);
};





int main(int argc, char *argv[]) {
    unsigned int M;
    unsigned int N;

    cout << "input log_2(number of grid points)" << endl;
    cin >> M;

    cout << "input test = 1,2, other" << endl;
    int test;
    cin >> test;

    cout << "input particle refinement factor" << endl;
    unsigned int cfactor;
    cin >> cfactor;

    cout << "enter stopping time" << endl;
    double timeStop;
    cin >> timeStop;


    ParticleSet p;
    N = Power(2, M);
    double h = 1. / N;
    double hp = h / cfactor; // pow(h,4./3.);
    int Np = 1. / hp;
    hp = 1. / Np;


    double delta = h;
    int pcfactor = 4 / cfactor;
    if (pcfactor < 1) pcfactor = 1;

    cout << "number of particles per cell = " << h * h / hp / hp << endl;


    shared_ptr<CutoffKernel> cutkptr = shared_ptr<CutoffKernel>(new CutoffKernel(h, delta));
    shared_ptr<ConvKernel> convkerptr = dynamic_pointer_cast<ConvKernel>(cutkptr);
    p.m_hockney.define(convkerptr, h, M);
    p.m_dx = h;


    array<double, DIM> lowCorner;
    if (test == 1) {

        // M = 6
        p.resize(1);
        Particle &particle = p.getParticle(0);
        particle.m_x[0] = .49;
        particle.m_x[1] = .24;
        particle.strength = 1. / h / h;

    } else if (test == 2) {

        // M = 6
        p.resize(2);

        Particle &particle = p.getParticle(0);
        particle.m_x[0] = .5;
        particle.m_x[1] = .5;
        particle.strength = 1. / h / h;

        Particle &particle2 = p.getParticle(1);
        particle2.m_x[0] = .5;
        particle2.m_x[1] = .25;
        particle2.strength = 0;

    } else if (test == 3) {
        // M = 6
        p.resize(2);

        Particle &particle = p.getParticle(0);
        particle.m_x[0] = .5;
        particle.m_x[1] = .75;
        particle.strength = 1. / h / h;

        Particle &particle2 = p.getParticle(1);
        particle2.m_x[0] = .5;
        particle2.m_x[1] = .25;
        // Opposite-sign vortices move together in a straight line (translational motion)
        particle2.strength = 1. / h / h;

    } else if (test == 4) {
        // Two-patch problem. Vertical
        // M = 7
        // Two closely spaced vortex patches (collections of vortex particles) with the same vorticity sign
        // may merge into a single vortex due to mutual attraction and diffusion. This is a classic test for
        // vortex methods.
        array<double, DIM> xp;
        assert(Np == 256 && cfactor == 2 && N == 128 && hp == 1. / 256 && h == 1. / 128 &&
               "Np = 256, cfactor = 2; N = 128, Np = 1. / hp = cfactor / h = cfactor * N;");
        for (int i = 0; i < Np; i++) {
            xp[0] = i * hp;
            for (int j = 0; j < Np; j++) {
                xp[1] = j * hp;
                double dist1 = sqrt(pow(xp[0] - .375, 2) + pow(xp[1] - .5, 2));
                double dist2 = sqrt(pow(xp[0] - .625, 2) + pow(xp[1] - .5, 2));
                if ((dist1 < .12) | (dist2 < .12)) {
                    Particle part;
                    part.m_x[0] = xp[0];
                    part.m_x[1] = xp[1];
                    part.strength = hp * hp / h / h;
                    p.addParticle(part);
                }
            }
        }

    } else if (test == 5) {
        // Two-patch problem. Horizontal
        // M = 7
        array<double, DIM> xp;
        assert(Np == 256 && cfactor == 2 && N == 128 && hp == 1. / 256 && h == 1. / 128 &&
               "Np = 256, cfactor = 2; N = 128, Np = 1. / hp = cfactor / h = cfactor * N;");
        for (int i = 0; i < Np; i++) {
            xp[0] = i * hp;
            for (int j = 0; j < Np; j++) {
                xp[1] = j * hp;
                double dist1 = sqrt(pow(xp[0] - .5, 2) + pow(xp[1] - .375, 2));
                double dist2 = sqrt(pow(xp[0] - .5, 2) + pow(xp[1] - .625, 2));
                if ((dist1 < .12) | (dist2 < .12)) {
                    Particle part;
                    part.m_x[0] = xp[0];
                    part.m_x[1] = xp[1];
                    part.strength = hp * hp / h / h;
                    p.addParticle(part);
                }
            }
        }

    } else if (test == 6) {
        // Two-patch problem. Less dense grid with longer time window
        // M = 10
        array<double, DIM> xp;
        Np = 256;
        assert(Np == 256 && cfactor == 2 && "Np = 256, cfactor = 2;");
        for (int i = 0; i < Np; i++) {
            xp[0] = i * hp;
            for (int j = 0; j < Np; j++) {
                xp[1] = j * hp;
                double dist1 = sqrt(pow(xp[0] - .5, 2) + pow(xp[1] - .375, 2));
                double dist2 = sqrt(pow(xp[0] - .5, 2) + pow(xp[1] - .625, 2));
                if ((dist1 < .12) | (dist2 < .12)) {
                    Particle part;
                    part.m_x[0] = xp[0];
                    part.m_x[1] = xp[1];
                    part.strength = hp * hp / h / h;
                    p.addParticle(part);
                }
            }
        }

    } else if (test == 7) { // M = 7
        // Vortex Sheet Roll-Up (Kelvin-Helmholtz Instability)
        // @see More Kelvin-Helmholtz (The Experiment) - Sixty Symbols
        // A vortex sheet (a line of closely spaced vortices) is unstable and rolls up into larger vortices
        // due to the Kelvin-Helmholtz instability.
        Np = N / cfactor;
        double x_start = 0.05; // Place particles in x in [0.05, 0.95] to avoid boundaries
        double x_end = 0.95;
        double width = x_end - x_start;
        hp = width / Np;
        int layers = 9;
        p.resize(Np * layers);

        // The vorticity field resembles a thin, wavy sheet centered around y = 0.5
        // with a sinusoidal perturbation in the y-direction.
        // As time evolves, the perturbation grows and the sheet rolls up into larger vortices.
        // Each vortex is a tightly wound spiral of particles, with high vorticity concentrated at the centers
        // and filamentary structures connecting them.

        double perturbation = 0.01; // Perturbation amplitude
        for (int i = 0; i < Np; i++) {
            for (int j = 0; j < layers; j++) {
                Particle &part = p.getParticle(i + j * Np);
                part.m_x[0] = x_start + i * hp; // x from 0.05 to < 0.95
                part.m_x[1] = 0.45 + perturbation * sin(2 * M_PI * (part.m_x[0] - x_start) / width);
                part.strength = hp / h / h / 5; // Uniform vorticity
                part.m_x[1] += j * 0.0008;      // Offset in y for each layer
            }
        }
    } else if (test == 8) {
        // @see Leapfrogging Vortex Rings - Reformulated Vortex Particle Method
        // @see The Science of Vortex Rings

        // Two vortex pairs (each with opposite-signed vortices) are initialized such that one pair passes
        // through the other, repeating periodically in a leapfrogging motion.

        p.resize(4);
        // Pair 1: Vortices at (0.75, 0.55) and (0.75, 0.45)
        Particle &p1 = p.getParticle(0);
        p1.m_x[0] = 0.75;
        p1.m_x[1] = 0.55;
        p1.strength = .1 / h / h;

        Particle &p2 = p.getParticle(1);
        p2.m_x[0] = 0.75;
        p2.m_x[1] = 0.45;
        p2.strength = -0.1 / h / h;

        // Pair 2: Vortices at (0.85, 0.55) and (0.85, 0.45)
        Particle &p3 = p.getParticle(2);
        p3.m_x[0] = 0.85;
        p3.m_x[1] = 0.55;
        p3.strength = 0.1 / h / h;

        Particle &p4 = p.getParticle(3);
        p4.m_x[0] = 0.85;
        p4.m_x[1] = 0.45;
        p4.strength = -.1 / h / h;

    } else if (test == 9) {

        // Leapfrogging Vortex Rings - Two pairs of vortex patches
        int particles_per_patch = 9;        ///< 3x3 grid per vortex
        p.resize(4 * particles_per_patch);  ///< 4 vortices * 9 particles
        double patch_radius = 0.01;         ///< Patch size
        double total_strength = .1 / h / h; ///< Strength per vortex
        double strength_per_particle = total_strength / particles_per_patch;

        // Grid for patch particles (3x3)
        double offsets[3] = {-patch_radius, 0, patch_radius};

        int idx = 0;
        for (int v = 0; v < 2; v++) { // Patches centered at (x1, y2) and (x1, y1)
            double x_center = 0.7;
            double y_center = 0.55 - v * 0.2;
            double sign = (v == 0) ? 1.0 : -1.0;
            for (int i = 0; i < 3; i++) {
                for (int j = 0; j < 3; j++) {
                    Particle &part = p.getParticle(idx);
                    part.m_x[0] = x_center + offsets[i];
                    part.m_x[1] = y_center + offsets[j];
                    part.strength = sign * strength_per_particle;
                    idx++;
                }
            }
        }

        for (int v = 0; v < 2; v++) { // Patches centered at (x2, y2) and (x2, y1)
            double x_center = 0.8;
            double y_center = 0.55 - v * 0.2;
            double sign = (v == 0) ? 1.0 : -1.0;
            for (int i = 0; i < 3; i++) {
                for (int j = 0; j < 3; j++) {
                    Particle &part = p.getParticle(idx);
                    part.m_x[0] = x_center + offsets[i];
                    part.m_x[1] = y_center + offsets[j];
                    part.strength = sign * strength_per_particle;
                    idx++;
                }
            }
        }

    } else if (test == 10) {

        // Periodic Vortex Array (Karman Vortex Street Approximation)
        // @see Von Karman Street of Vortices (2) Experimental observations

        int numVortices = 64; // Two rows of 4 vortices
        int halfNumVortices = numVortices / 2;
        p.resize(numVortices);

        double dx = 0.01; // Horizontal spacing for stability
        double dy = 0.03;

        double x_start = 0.5;
        double y_top = 0.55;
        double y_bottom = y_top - dy;
        double strength = 0.006 / h / h;

        for (int i = 0; i < halfNumVortices; i++) {
            // Top row: Positive vortices
            Particle &p1 = p.getParticle(i);
            p1.m_x[0] = x_start + i * dx;
            p1.m_x[1] = y_top;
            p1.strength = strength;

            // Bottom row: Negative vortices, staggered
            Particle &p2 = p.getParticle(i + halfNumVortices);
            p2.m_x[0] = x_start + (i + 0.5) * dx;
            p2.m_x[1] = y_bottom;
            p2.strength = -strength;
        }

    } else {
        array<double, DIM> xp;

        for (int i = 0; i < Np; i++) {
            xp[0] = i * hp;

            for (int j = 0; j < Np; j++) {
                xp[1] = j * hp;
                double dist1 = sqrt(pow(xp[0] - .375, 2) + pow(xp[1] - .5, 2));
                double dist2 = sqrt(pow(xp[0] - .625, 2) + pow(xp[1] - .5, 2));
                if ((dist1 < .12) | (dist2 < .12)) {
                    Particle part;
                    part.m_x[0] = xp[0];
                    part.m_x[1] = xp[1];
                    part.strength = hp * hp / h / h;
                    p.addParticle(part);
                }
            }
        }
    }

    double dx = 1. / N;
    cout << "number of particles = " << p.getSize() << endl;

    Point low = getZeros();
    Point high = getOnes() * N;
    DBox bx(low, high);
    p.m_box = bx;

    lowCorner[0] = 0.;
    lowCorner[1] = 0.;

    ParticleShift kIn, kOut;
    kIn.init(p);
    kOut.init(p);
    kIn.setToZero();

    ParticleVelocities pv();
    double time = 0.;
    double dt = 140 * .025 / N;
    int m = 500000;

    RK4<ParticleSet, ParticleVelocities, ParticleShift> integrator;
    outField(p, pcfactor);
    PWrite(&p);



    for (int i = 0; i < m; i++) {

        integrator.advance(time, dt, p);

        time = time + dt;

#define TEST89
#ifdef TEST89
        for (int k = 0; k < p.getSize(); k++) {
            Particle &part = p.getParticle(k);
            part.m_x[0] = std::max(0.0, std::min(1.0, part.m_x[0]));
            part.m_x[1] = std::max(0.0, std::min(1.0, part.m_x[1]));
        }
#endif

        if (i % 10 == 0) {
            cout << "time = " << time << "  dt " << dt << endl;
        }

        outField(p, pcfactor);
        PWrite(&p);


        if (time >= timeStop) {
            cout << "Early stop at time = " << time << endl;
            break;
        }
    }
}
