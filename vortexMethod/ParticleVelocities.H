#ifndef _PARTICLEVELOCITIES_H_
#define _PARTICLEVELOCITIES_H_

#include "DBox.H"
#include "Hockney.H"
#include "ParticleSet.H"
#include "RectMDArray.H"

#include <cassert>
#include <cmath>
#include <iostream>
#include <vector>

using namespace std;

class ParticleVelocities {
  public:
    /**
     * @brief Default constructor
     */
    ParticleVelocities();

    /**
     * @brief () function: implementing one piece of the Runge-Kutta 4 (RK4) integration method.
     *  This function computes the velocity of particles based on their current state and time, which will be
     *  used to update their positions over time.
     *
     *  The stages of RK4 all require the calculation of quantities of the form k := dt * F(t + dt, X + k)
     *
     *  In this function:
     *  1) Deposit particle vorticity onto a grid (taking into account the displacement in a_k)
     *  2) Solve the Poisson equation using the Hockney method to get the stream function
     *  3) Compute the velocity field on the grid using finite differences
     *  4) Interpolate the velocities back to particles
     *  5) Scale the velocities by dt and store them in a_k
     *
     * @param a_k The ParticleShift object representing the new state for k, used both
     *  to store the result dt·U and as the displacement to apply to the particles before evaluating.
     *  This displacement helps RK4 evaluate the system at intermediate states.
     *
     * @param a_time The current time to evaluate the function t + dt
     * @param dt The time step size
     * @param a_state The ParticleSet object representing the current state of the particles
     */
    void operator()(ParticleShift &a_k, const double &a_time, const double &dt, ParticleSet &a_state);
    // This member function computes the k's induced on a set of particles
    // by all of the particles in the input ParticleSet
    // displaced by the input k.


    /**
     * @brief Compute the displacement increment for each particle in the ParticleSet
     *  With a modularized version of the workflow
     */
    void computeDisplacementIncrement_detailed_v0(ParticleShift &a_k, const double &a_time, const double &dt,
                                                  ParticleSet &a_state);
    /**
     * @brief Compute the displacement increment for each particle in the ParticleSet
     */
    void computeDisplacementIncrement(ParticleShift &a_k, const double &a_time, const double &dt,
                                      ParticleSet &a_state);
};
#endif
