#ifndef _PARTICLEVELOCITIES_H_
#define _PARTICLEVELOCITIES_H_
#include <iostream>
#include <cassert>
#include <cmath>
#include <vector>
#include "DBox.H"
#include "RectMDArray.H"
#include "Hockney.H"
#include "ParticleSet.H"
using namespace std;
class ParticleVelocities
{
public:
  ParticleVelocities();
  void operator()(ParticleShift& a_k, 
                     const double& a_time, const double& dt, 
                     ParticleSet& a_state);
};
#endif
