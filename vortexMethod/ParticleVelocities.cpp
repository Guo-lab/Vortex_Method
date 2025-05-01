#include "ParticleVelocities.H"

#include "DBox.H"
#include "Point.H"
#include "RectMDArray.H"




ParticleVelocities::ParticleVelocities() {}




/**
 * @brief Doing the particle velocity computation
 *
 * @note
 * (1 − s0k) (1 − s1k): Fraction of the cell diagonally opposite the particle
 * (closer to the lower-left corner).
 *      s0k (1 − s1k) : Fraction closer to the lower-right corner.
 *      (1 − s0k) s1k : Fraction closer to the upper-left corner.
 *            s0k s1k : Fraction closer to the upper-right corner.
 *
 * s0k is the fractional distance in the x-direction from the lower-left.
 * s1k is the fractional distance in the y-direction from the lower-left.
 *
 * @details
 * @li Step 1: Deposit vorticity to grid
 *      - Initialize the vorticity grid to zero
 *      - For each particle, compute the displaced position
 *      - Find the nearest grid point below the particle and determine the cell containing the particle
 * @li Step 2: Solve Poisson equation
 * @li Step 3: Compute velocity field on grid
 *      - For each grid point, compute the velocity using central differences
 *      - Store the velocity in the U_grid
 * @li Step 4: Interpolate velocities to particles
 *      - For each particle, find the surrounding grid points
 *      - Compute the velocity at the particle's position using bilinear interpolation
 *      - Update the ParticleShift with the computed velocity
 */
void ParticleVelocities::computeDisplacementIncrement(ParticleShift &a_k, const double &a_time,
                                                      const double &dt, ParticleSet &a_state) {
    // Grid parameters
    double h = a_state.getH();
    RectMDArray<double> omega_grid(a_state.m_box);         ///< Vorticity grid
    RectMDArray<double> psi_grid(a_state.m_box);           ///< Stream function grid
    RectMDArray<array<double, DIM>> U_grid(a_state.m_box); ///< Velocity grid

    // Const variables
    DBox domain = omega_grid.getDBox();       ///< Domain of the grid
    Point e0 = getUnitv(0), e1 = getUnitv(1); ///< Unit vectors in the grid
    array<double, DIM> x_displaced;           ///< Displaced position of the particles
    std::array<int, DIM> i_k;                 ///< Grid around the particles
    std::array<double, DIM> s_k;              ///< Distance from the particles to the surrounding grid points

    // <i> Deposit vorticity to grid
    omega_grid.setVal(0.0);
    psi_grid.setVal(0.0);
    for (size_t k = 0; k < a_state.getSize(); k++) {
        for (int d = 0; d < DIM; d++) x_displaced[d] = a_state.getParticle(k).m_x[d] + a_k.getShift(k).m_x[d];
        std::transform(x_displaced.begin(), x_displaced.end(), i_k.begin(),
                       [h](double x) { return static_cast<int>(std::floor(x / h)); });
        std::transform(x_displaced.begin(), x_displaced.end(), i_k.begin(), s_k.begin(),
                       [h](double x, int i) { return (x - i * h) / h; });
        Point p_ik(i_k);
        double strength = a_state.getParticle(k).strength;
        if (domain.contains(p_ik)) omega_grid[p_ik] += strength * (1.0 - s_k[0]) * (1.0 - s_k[1]);
        if (domain.contains(p_ik + e0)) omega_grid[p_ik + e0] += strength * s_k[0] * (1.0 - s_k[1]);
        if (domain.contains(p_ik + e1)) omega_grid[p_ik + e1] += strength * (1.0 - s_k[0]) * s_k[1];
        if (domain.contains(p_ik + e0 + e1)) omega_grid[p_ik + e0 + e1] += strength * s_k[0] * s_k[1];
    }

#ifdef DEBUG_LOG
    std::cout << "Deposited vorticity onto grid." << std::endl;
#endif

    // <ii> Solve Poisson equation
    a_state.m_hockney.convolve(omega_grid);
    psi_grid = omega_grid;

    // <iii> Compute velocity field on grid
    array<double, DIM> zeros = {0};
    U_grid.setVal(zeros);
    DBox interiorDomain = domain.grow(-1);
    for (Point pt = interiorDomain.getLowCorner(); interiorDomain.notDone(pt); interiorDomain.increment(pt)) {
        U_grid[pt][0] = (psi_grid[pt + getUnitv(1)] - psi_grid[pt - getUnitv(1)]) / (2.0 * h);  // u = dψ/dy
        U_grid[pt][1] = -(psi_grid[pt + getUnitv(0)] - psi_grid[pt - getUnitv(0)]) / (2.0 * h); // v = -dψ/dx
    }

#ifdef DEBUG_LOG
    std::cout << "Velocity field computed." << std::endl;
#endif

    // <iv> Interpolate velocities to particles and update ParticleShift (DIM 0 is x, DIM 1 is y)
    for (int k = 0; k < a_state.getSize(); k++) {
        for (int d = 0; d < DIM; d++) x_displaced[d] = a_state.getParticle(k).m_x[d] + a_k.getShift(k).m_x[d];
        std::transform(x_displaced.begin(), x_displaced.end(), i_k.begin(),
                       [h](double x) { return static_cast<int>(std::floor(x / h)); });
        std::transform(x_displaced.begin(), x_displaced.end(), i_k.begin(), s_k.begin(),
                       [h](double x, int i) { return (x - i * h) / h; });
        array<double, DIM> U_k = {0.0};
        Point p_ik(i_k);

        if (domain.contains(p_ik) && domain.contains(p_ik + e0) && domain.contains(p_ik + e1) &&
            domain.contains(p_ik + e0 + e1)) {
            Point p_lr = p_ik + e0, p_ul = p_ik + e1, p_ur = p_ik + e0 + e1;
            for (int d = 0; d < DIM; d++) { // U_k = (u, v) = (dψ/dy, -dψ/dx)
                U_k[d] = U_grid[p_ik][d] * (1 - s_k[0]) * (1 - s_k[1]) +
                         U_grid[p_lr][d] * s_k[0] * (1 - s_k[1]) + U_grid[p_ul][d] * (1 - s_k[0]) * s_k[1] +
                         U_grid[p_ur][d] * s_k[0] * s_k[1];
            }
        } else {
            // Handle boundary case (e.g., set velocity to zero or log a warning)
            U_k[0] = 0.0;
            U_k[1] = 0.0;
        }
        for (int d = 0; d < DIM; d++) a_k.getShift(k).m_x[d] = dt * U_k[d];
#ifdef DEBUG_LOG
        std::cout << "PV.cpp Particle " << k << ": Shift = (" << shift.m_x[0] << ", " << shift.m_x[1] << ")"
                  << std::endl;
#endif
    }

#ifdef DEBUG_LOG
    std::cout << "Interpolated velocities to particles." << std::endl;
#endif
}



// Helper functions
void depositParticlesToGrid(RectMDArray<double> &a_gridVorticity, const ParticleSet &a_state,
                            const ParticleShift &a_k, const double &dx) {
    for (size_t k = 0; k < a_state.getSize(); k++) {
        array<double, DIM> x_shifted;
        for (int d = 0; d < DIM; d++) {
            x_shifted[d] = a_state.getParticle(k).m_x[d] + a_k.getShift(k).m_x[d];
        }

        // Find nearest grid point below particle. Determine the cell containing the particle
        array<int, DIM> ik;
        array<double, DIM> sk;
        for (int d = 0; d < DIM; d++) {
            ik[d] = floor(x_shifted[d] / dx);
            sk[d] = (x_shifted[d] - ik[d] * dx) / dx;
        }

        // Distribute vorticity to the four nearest grid points in 2D case
        double weight00 = (1.0 - sk[0]) * (1.0 - sk[1]);
        double weight10 = sk[0] * (1.0 - sk[1]);
        double weight01 = (1.0 - sk[0]) * sk[1];
        double weight11 = sk[0] * sk[1];

        // Add contributions to grid vorticity
        Point p_ik(ik), e0 = getUnitv(0), e1 = getUnitv(1);
        double particle_strength = a_state.getParticle(k).strength;

        DBox domain = a_gridVorticity.getDBox();
        if (domain.contains(p_ik)) a_gridVorticity[p_ik] += particle_strength * weight00;
        if (domain.contains(p_ik + e0)) a_gridVorticity[p_ik + e0] += particle_strength * weight10;
        if (domain.contains(p_ik + e1)) a_gridVorticity[p_ik + e1] += particle_strength * weight01;
        if (domain.contains(p_ik + e0 + e1)) a_gridVorticity[p_ik + e0 + e1] += particle_strength * weight11;
    }
}



void computeVelocityField(RectMDArray<array<double, DIM>> &a_velocityField,
                          const RectMDArray<double> &a_streamFunc, const DBox &a_domain, const double &dx) {
    // Create interior domain (shrinked by 1)
    // excluding boundary points for central differences
    DBox interiorDomain = a_domain.grow(-1);

    for (Point pt = interiorDomain.getLowCorner(); interiorDomain.notDone(pt); interiorDomain.increment(pt)) {
        Point ip0 = pt + getUnitv(0);
        Point im0 = pt - getUnitv(0);
        Point ip1 = pt + getUnitv(1);
        Point im1 = pt - getUnitv(1);
        // Compute velocity using central differences U = (dψ/dy, -dψ/dx)
        a_velocityField[pt][0] = (a_streamFunc[ip1] - a_streamFunc[im1]) / (2.0 * dx);
        a_velocityField[pt][1] = -(a_streamFunc[ip0] - a_streamFunc[im0]) / (2.0 * dx);
    }
}


void interpolateVelocitiesToParticles(ParticleShift &a_k,
                                      const RectMDArray<array<double, DIM>> &a_velocityField,
                                      const ParticleSet &a_state, const double &dx, const double &dt) {
    for (size_t k = 0; k < a_state.getSize(); k++) {
        array<double, DIM> x_shifted;
        for (int d = 0; d < DIM; d++) {
            x_shifted[d] = a_state.getParticle(k).m_x[d] + a_k.getShift(k).m_x[d];
        }

        array<int, DIM> ik;
        array<double, DIM> sk;
        for (int d = 0; d < DIM; d++) {
            ik[d] = floor(x_shifted[d] / dx);
            sk[d] = (x_shifted[d] - ik[d] * dx) / dx;
        }

        double weight00 = (1.0 - sk[0]) * (1.0 - sk[1]);
        double weight10 = sk[0] * (1.0 - sk[1]);
        double weight01 = (1.0 - sk[0]) * sk[1];
        double weight11 = sk[0] * sk[1];

        Point p_ik(ik), e0 = getUnitv(0), e1 = getUnitv(1);

        // Interpolate velocity to particle position
        array<double, DIM> velocity = {0.0, 0.0};
        for (int d = 0; d < DIM; d++) {
            // Compute interpolated velocity
            DBox domain = a_velocityField.getDBox();
            if (domain.contains(p_ik)) velocity[d] += weight00 * a_velocityField[p_ik][d];
            if (domain.contains(p_ik + e0)) velocity[d] += weight10 * a_velocityField[p_ik + e0][d];
            if (domain.contains(p_ik + e1)) velocity[d] += weight01 * a_velocityField[p_ik + e1][d];
            if (domain.contains(p_ik + e0 + e1)) velocity[d] += weight11 * a_velocityField[p_ik + e0 + e1][d];
        }

        // Update a_k with dt * velocity
        for (int d = 0; d < DIM; d++) {
            a_k.getShift(k).m_x[d] = dt * velocity[d];
        }
    }
}

void ParticleVelocities::computeDisplacementIncrement_detailed_v0(ParticleShift &a_k, const double &a_time,
                                                                  const double &dt, ParticleSet &a_state) {
    DBox domain = a_state.m_box;
    double dx = a_state.m_dx;

    RectMDArray<double> vorticityGrid(domain);
    RectMDArray<array<double, DIM>> velocityField(domain);
    vorticityGrid.setVal(0.0);
    array<double, DIM> zeros = {0};
    velocityField.setVal(zeros);

    // Deposit the vorticity from particles (with shifts) onto the grid
    depositParticlesToGrid(vorticityGrid, a_state, a_k, dx);

    // Solve the Poisson equation using Hockney algorithm
    a_state.m_hockney.convolve(vorticityGrid);

    // Compute the velocity field on the grid using finite differences
    computeVelocityField(velocityField, vorticityGrid, domain, dx);

    // Interpolate velocities back to particles and scale by dt
    interpolateVelocitiesToParticles(a_k, velocityField, a_state, dx, dt);
}





void ParticleVelocities::operator()(ParticleShift &a_k, const double &a_time, const double &dt,
                                    ParticleSet &a_state) {
    // #define DEBUG_WORKFLOW

#ifdef DEBUG_WORKFLOW
    computeDisplacementIncrement_detailed_v0(a_k, a_time, dt, a_state);
#else
    computeDisplacementIncrement(a_k, a_time, dt, a_state);
#endif
}
