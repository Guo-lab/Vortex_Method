
#ifndef _RK4_H_
#define _RK4_H_

#include <iostream>

///  generic explicit RK4 algorithm
/**  class X is the type of your state data. class F provides an evaluation
     of the X' for a specified state of X and the time.  If your X' is
     determined entirely from the current state values then F can ignore.
     F implements the function:
     void operator(X& k, double a_time, double a_dt, const X& a_state, const X& a_shift);

     F interprets the real state as a_state+a_shift

     Note that F includes the customary factor of dT in the RK4 algorithm.
     (i.e. this algorithm assumes k_i = F(...) NOT k_i = dt*F(...)

     like RK4<ParticleSet, ParticleVelocities, ParticleShift> integrator;
*/
template <class X, class F, class dX> class RK4 {
  public:
    void advance(double a_time, double a_dt, X &a_state);

  protected:
    dX m_k;
    dX m_delta;
    F m_f;
};

//==============================================

template <class X, class F, class dX> void RK4<X, F, dX>::advance(double a_time, double a_dt, X &a_state) {
    // this looks odd, but it is a common programming error you can't write "1/6* m_k";
    // and get what you want. the first term is evaluated as integer division, resulting in a zero.

    /**
     * @ref https://lpsa.swarthmore.edu/NumInt/NumIntFourth.html
     * @brief RK4 algorithm
     *  m       := 1/6 k1 + 1/3 k2 + 1/3 k3 + 1/6 k4
     *  a_state := current estimate of k
     *  m_k     := new estimate of k
     */

#ifdef DEBUG_RK4
    for (int k = 0; k < a_state.getSize(); k++) {
        std::cout << "Begin RK4 Particle " << k << ": (" << a_state.getParticle(k).m_x[0] << ", "
                  << a_state.getParticle(k).m_x[1] << ")" << std::endl;
    }
#endif
    double sixth = 1, third = 1, half = 1;
    sixth /= 6;
    third /= 3;
    half /= 2;

    m_delta.init(a_state);
    m_k.init(a_state);

#ifdef DEBUG_RK4
    std::cout << "m_delta and m_k initialized" << std::endl;
#endif
    m_f(m_k, a_time, a_dt, a_state); // compute k1 = f(tn, xn)
#ifdef DEBUG_RK4
    std::cout << "k1 done" << std::endl;
    for (int k = 0; k < a_state.getSize(); k++) {
        std::cout << "RK4 Particle " << k << ": Shift = (" << m_k.getShift(k).m_x[0] << ", "
                  << m_k.getShift(k).m_x[1] << ")" << std::endl;
    }
#endif

    m_delta.increment(sixth, m_k); // add dt k1/6
    m_k *= half;                   // shift: x + dt k1/2

    m_f(m_k, a_time + half * a_dt, a_dt, a_state); // compute k2 = f(tn + dt/2, xn + k1*dt/2)
#ifdef DEBUG_RK4
    std::cout << "k2 done" << std::endl;
    for (int k = 0; k < a_state.getSize(); k++) {
        std::cout << "RK4 Particle " << k << ": Shift = (" << m_k.getShift(k).m_x[0] << ", "
                  << m_k.getShift(k).m_x[1] << ")" << std::endl;
    }
#endif


    m_delta.increment(third, m_k); // add dt k2/3
    m_k *= half;                   // shift: x + dt k2/2

    m_f(m_k, a_time + half * a_dt, a_dt, a_state); // conpute k3 = f(tn + dt/2, xn + k2*dt/2)
#ifdef DEBUG_RK4
    std::cout << "k3 done" << std::endl;
    for (int k = 0; k < a_state.getSize(); k++) {
        std::cout << "RK4 Particle " << k << ": Shift = (" << m_k.getShift(k).m_x[0] << ", "
                  << m_k.getShift(k).m_x[1] << ")" << std::endl;
    }
#endif


    m_delta.increment(third, m_k); // add dt k3/3

    m_f(m_k, a_time + a_dt, a_dt, a_state); // compute k4 = f(tn + dt, xn + k3*dt)
#ifdef DEBUG_RK4
    std::cout << "k4 done" << std::endl;
    for (int k = 0; k < a_state.getSize(); k++) {
        std::cout << "RK4 Particle " << k << ": Shift = (" << m_k.getShift(k).m_x[0] << ", "
                  << m_k.getShift(k).m_x[1] << ")" << std::endl;
    }
#endif

    m_delta.increment(sixth, m_k); // add dt k4/6

    // y*(t0) + (k1/6 + k2/3 + k3/3 + k4/6) dt.
    a_state.increment(m_delta);
#ifdef DEBUG_RK4
    for (int k = 0; k < a_state.getSize(); k++) {
        std::cout << "End   RK4 Particle " << k << ": (" << a_state.getParticle(k).m_x[0] << ", "
                  << a_state.getParticle(k).m_x[1] << ")" << std::endl;
    }
#endif
}

#endif
