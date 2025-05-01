#ifndef _PARTICLESET_H_
#define _PARTICLESET_H_

#include <vector>
#include "Particle.H"
#include "RectMDArray.H"
#include "Hockney.H"

using namespace std;

class ParticleSet;
class ParticleShift
{
public:
  vector<DX> m_particles;
  void init(const ParticleSet& a_particles);
  void increment(double a_scale, const ParticleShift& a_rhs);
  void operator*=(double a_scale);
  void setToZero();
};
class ParticleSet
{
public:
  //ParticleSet(Box& a_box,double& a_dx,double a_lowCorner[DIM],int a_M);
  ParticleSet(shared_ptr<ConvKernel>& a_kerptr,Box& a_box,double& a_dx, 
              array<Real, DIM>& a_lowCorner,int a_M);
  ParticleSet(){};
  vector<Particle> m_particles;
  double m_dx;
  Box m_box;
  double m_lowCorner[DIM];
  Hockney m_hockney;

  void increment(const ParticleShift& a_shift);
};


#endif
