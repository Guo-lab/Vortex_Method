#include "ParticleSet.H"

/// @brief Initializes the ParticleShift object with the given ParticleSet
/// @param a_particles
// void ParticleShift::init(const ParticleSet &a_particles) {
//     size_t n_particles = a_particles.getSize();
//     m_particles.resize(n_particles);
//     for (unsigned int i = 0; i < (unsigned int)n_particles; i++) {
//         std::copy(a_particles.getParticle(i).m_x.begin(), a_particles.getParticle(i).m_x.end(),
//                   m_particles[i].m_x.begin());
//     }
// }
void ParticleShift::init(const ParticleSet &a_particles) {
    size_t n_particles = a_particles.getSize();
    m_particles.resize(n_particles);
    for (unsigned int i = 0; i < (unsigned int)n_particles; i++) {
        // Initialize shifts to zero, not particle positions
        m_particles[i].m_x.fill(0.0);
    }
}


/// @brief Initializes the ParticleSet object with the given ConvKernel, DBox, dx, low corner, and M
/// @param a_kerptr
/// @param a_box
/// @param a_dx
/// @param a_lowCorner
/// @param a_M
ParticleSet::ParticleSet(std::shared_ptr<ConvKernel> &a_kerptr, DBox &a_box, double &a_dx,
                         std::array<double, DIM> &a_lowCorner, int a_M) {
    m_hockney.define(a_kerptr, a_dx, a_M);
    m_dx = a_dx;
    m_box = a_box;
    for (unsigned int i = 0; i < DIM; i++) {
        m_lowCorner[i] = a_lowCorner[i];
    }
}
