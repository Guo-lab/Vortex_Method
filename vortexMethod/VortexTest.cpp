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
        particle2.strength = 1. / h / h;
    } else if (test == 4) {
        // Two-patch problem.
        // M = 7
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
        // Two-patch problem.
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
        // Two-patch problem.
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
