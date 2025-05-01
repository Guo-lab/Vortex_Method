#ifndef _PARTICLE_H_
#define _PARTICLE_H_
#include <cstdio>
#include <cmath>
#include <cassert>
#include <array>
using namespace std;
class DX
{
public:
  DX()
  {
    for(unsigned int i=0; i<DIM; i++) 
      {
	m_x[i]=0;
      }
  } 
  array<double, DIM> m_x;
  inline void increment(double a_scale, const DX& a_rhs)
  {
    for(unsigned int i=0; i<DIM; i++) 
      {
	m_x[i]+=a_scale*a_rhs.m_x[i];
      }
  }
  inline void operator*=(double a_scale)
  {
    for(unsigned int i=0; i<DIM; i++) 
      {
	m_x[i]*=a_scale;
      }
  }
};
class Particle
{
public:
  array<double, DIM> m_x;
  double strength;
  inline void increment(const DX& a_shift)
  {
    for(unsigned int i=0; i<DIM; i++) 
      {
	m_x[i]+=a_shift.m_x[i];
      }
  }
};

#endif
