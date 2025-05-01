#ifndef _PARTICLESET_H_
#define _PARTICLESET_H_

#include "Hockney.H"
#include "Particle.H"
#include "RectMDArray.H"

#include <vector>

using namespace std;

class ParticleSet;
class ParticleVelocities;

/**
 * @brief Represents shifts/displacements of particles.
 *  Used during integration to store intermediate derivatives.
 */
class ParticleShift {
  public:
    /**
     * @brief Default constructor
     *  Initializes the vector of particles to an empty state
     */
    ParticleShift() {};

    /**
     * @brief Initializes the ParticleShift object with the given ParticleSet
     *
     * @param a_particles The ParticleSet to initialize from
     */
    void init(const ParticleSet &a_particles);

    /**
     * @brief Increment this ParticleShift by a scaled version of another ParticleShift
     *  Performs the operation: this += scale * rhs
     *
     * @param a_scale The scaling factor to apply to rhs
     * @param a_rhs The ParticleShift to scale and add to this ParticleShift
     */
    void increment(double a_scale, const ParticleShift &a_rhs) {
        assert(m_particles.size() == a_rhs.getSize());
        for (unsigned int i = 0; i < m_particles.size(); i++) {
            m_particles[i].increment(a_scale, a_rhs.getShift(i));
        }
    }

    /**
     * @brief Scale this ParticleShift by a scalar value
     *  Performs the operation: this *= scale
     *
     * @param a_scale The scaling factor to apply to this ParticleShift
     */
    void operator*=(double a_scale) {
        for (unsigned int i = 0; i < m_particles.size(); i++) {
            m_particles[i] *= a_scale;
        }
    }

    /**
     * @brief Set all components of this ParticleShift to zero
     *
     * This is used to reset the ParticleShift before integration
     */
    void setToZero() {
        for (unsigned int i = 0; i < m_particles.size(); i++) {
            m_particles[i].m_x.fill(0.0);
        }
    }

    /**
     * @brief Get the shift/displacement vector of particle at index i
     *
     * @param i The index of the particle
     * @return The shift vector for the particle at index i
     */
    DX &getShift(size_t i) { return m_particles[i]; }

    const DX &getShift(size_t i) const { return m_particles[i]; }

    /**
     * @brief Get the number of particles in this ParticleShift
     *
     * @return The number of particles
     */
    size_t getSize() const { return m_particles.size(); }

  private:
    /**
     * @brief Array holding one displacement per particle
     */
    vector<DX> m_particles;
};

/**
 * @brief Represents a set of particles with associated grid and solver.
 *  Contains the particles, their properties, and the computational domain
 */
class ParticleSet {
  public:
    /**
     * @brief Constructor
     */
    ParticleSet() {};

    /**
     * @brief Constructor with initialization of the ParticleSet with the given parameters
     *
     * @param a_kerptr The convolution kernel to use
     * @param a_box The computational domain
     * @param a_dx The grid spacing
     * @param a_lowCorner The lower corner coordinates of the computational domain
     * @param a_M The number of grid points in each dimension
     */
    ParticleSet(std::shared_ptr<ConvKernel> &a_kerptr, DBox &a_box, double &a_dx,
                std::array<double, DIM> &a_lowCorner, int a_M);

    /**
     * Increment the particles in the ParticleSet by a given shift from the ParticleShift
     *
     * @param a_shift The shift vector to apply to the particles
     */
    void increment(const ParticleShift &a_shift) {
        assert(m_particles.size() == a_shift.getSize());
        for (unsigned int i = 0; i < m_particles.size(); i++) {
            m_particles[i].increment(a_shift.getShift(i));
        }
    }

    /**
     * @brief Get the number of particles in the ParticleSet
     *
     * @return The number of particles
     */
    size_t getSize() const { return m_particles.size(); }

    /**
     * @brief Get the particles in the ParticleSet
     *
     * @param i The index of the particle
     * @return A vector of particles
     */
    Particle &getParticle(size_t i) { return m_particles[i]; }

    const Particle &getParticle(size_t i) const { return m_particles[i]; }

    /**
     * @brief Resize the particle container
     * @param size New number of particles
     */
    void resize(size_t size) { m_particles.resize(size); }

    /**
     * @brief Reserve capacity for particles
     * @param capacity Number of particles to reserve space for
     */
    void reserve(size_t capacity) { m_particles.reserve(capacity); }

    /**
     * @brief Add a particle to the ParticleSet
     *
     * @param particle The particle to add
     */
    void addParticle(const Particle &particle) { m_particles.push_back(particle); }

    /**
     * @brief Add a particle to the ParticleSet with given coordinates and strength
     * @param x The coordinates of the particle
     * @param strength The strength of the particle
     */
    void addParticle(const array<double, DIM> &coord, double strength) {
        Particle particle;
        particle.m_x = coord;
        particle.strength = strength;
        m_particles.push_back(particle);
    }

    double getH() const {
        assert(m_hockney.isDefined() && m_hockney.getH() == m_dx);
        return m_hockney.getH();
    }

    /**
     * @brief The Grid spacing
     */
    double m_dx;

    /**
     * @brief The computational domain
     */
    DBox m_box;

    /**
     * @brief The lower corner coordinates of the computational domain
     */
    double m_lowCorner[DIM];

    /**
     * @brief The Hockney solver for the convolution (for Poisson equation)
     */
    Hockney m_hockney;

  private:
    /**
     * @brief Collection of particles in the simulation
     */
    vector<Particle> m_particles;

    friend class ParticleVelocities;
};

#endif
