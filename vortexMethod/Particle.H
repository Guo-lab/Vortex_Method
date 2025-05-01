#ifndef _PARTICLE_H_
#define _PARTICLE_H_

#include <array>
#include <cassert>
#include <cmath>
#include <cstdio>

using namespace std;

#ifndef DIM
#define DIM 2
#endif

/**
 * @brief DX class represents a differential operator in the context of particle simulation.
 *  It holds an array of doubles representing the displacement vector.
 */
class DX {
  public:
    /**
     * @brief Array holding the DIM-dimensional vector components
     */
    array<double, DIM> m_x;

    /**
     * @brief Default constructor
     *
     * Initializes all components of the vector to zero
     */
    DX() : m_x() {}

    /**
     * @brief Increment this vector by a scaled version of another vector
     *
     * Performs the operation: this += scale * rhs
     *
     * @param a_scale The scaling factor to apply to rhs
     * @param a_rhs The vector to scale and add to this vector
     */
    inline void increment(double a_scale, const DX &a_rhs) {
        for (unsigned int i = 0; i < DIM; i++) {
            m_x[i] += a_scale * a_rhs.m_x[i];
        }
    }

    /**
     * @brief Scale this vector by a scalar value
     *
     * Performs the operation: this *= scale
     *
     * @param a_scale The scaling factor to apply to this vector
     */
    inline void operator*=(double a_scale) {
        for (unsigned int i = 0; i < DIM; i++) {
            m_x[i] *= a_scale;
        }
    }
};

/**
 * @brief Particle class represents a physical particle in the vortex simulation.
 *  It holds the position of the particle and its strength.
 */
class Particle {
  public:
    /**
     * @brief Array holding the DIM-dimensional position of the particle
     */
    array<double, DIM> m_x;

    /**
     * @brief Strength of the particle
     *
     * Represents the intensity/weights or influence of the particle in the simulation
     */
    double strength;

    /**
     * @brief Updates the particle's position by adding the components of the shift vector.
     *
     * @param a_shift The displacement vector to add to this particle's position
     */
    inline void increment(const DX &a_shift) {
        for (unsigned int i = 0; i < DIM; i++) {
            m_x[i] += a_shift.m_x[i];
        }
    }
};

#endif
