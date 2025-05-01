#include <cmath>
#include <complex>
#include <vector>
#include <cstdio>
#include <iostream>
using namespace std;
#include "PowerItoI.H"
#include "FFT1D.H"
#include "ConvKernel.H"
#include "CutoffKernel.H"
#include "FFTW1D.H"
#include "RectMDArray.H"
#include "WriteRectMDArray.H"
#include "Box.H"
#include "FFTMD.H"
#include "Hockney.H"
void hockneyTest(const int& a_M)
{
  int N = Power(2,a_M);
  int low[DIM],high[DIM],tuple[DIM];
  
  for (int dir = 0;dir <DIM;dir++)
    {
      low[dir] = 0;
      high[dir] = N;
    }
  Box b(low,high);
  for (int dir = 0;dir <DIM;dir++)
    {
      low[dir] = 0;
      high[dir] = 2*N;
    }
  Box bDouble(low,high);
  RectMDArray<double>  rho(b), rhoSave(b);
  rho.setVal(0.);
  rhoSave.setVal(0.);
  double h = 1./(N);
  rho.setVal(0.);
  double r0 = .25;
  int half[2] = {N/2,N/2};
  Point pth(half);
  //rho[pth] = 1./pow(h,DIM);
  for (Point pt = low; b.notDone(pt); b.increment(pt))
    {
      double r = sqrt(pow(pt[0]*h - .5,2) + pow(pt[1]*h - .5,2))/r0;
      if (r < 1.) 
        {
          rho[pt] = pow(cos(M_PI*r/2),6
);
          rhoSave[pt] = rho[pt];
        }
    }
  double delta = pow(h,.75);
  cout << "cutoff radius in cells = " << delta/h << endl;
  shared_ptr<CutoffKernel> cutkptr = 
    shared_ptr<CutoffKernel>(new CutoffKernel(h,delta));
  shared_ptr<ConvKernel> convkerptr = dynamic_pointer_cast<ConvKernel>(cutkptr);
  Hockney hockney(convkerptr,h,a_M);
  MDWrite("rho",rho);
  hockney.convolve(rho);
  MDWrite("convRho",rho);
  cout << "max value of G*rho " << rho[pth] << endl;
  Box bl = b.grow(-2);
  RectMDArray<double> LRho(bl);
  Point unit0 = getUnitv(0);
  Point unit1 = getUnitv(1);

  for (Point pt = bl.getLowCorner(); bl.notDone(pt); bl.increment(pt))
    {
      LRho[pt] = (-20*rho[pt] + 4*(rho[pt + unit0] + rho[pt-unit0] +
                                   rho[pt + unit1] + rho[pt-unit1]) +
                  rho[pt + unit0 + unit1] + rho[pt - unit0 + unit1] 
                  +rho[pt + unit0 - unit1] + rho[pt - unit0 - unit1])/6/h/h -
        (rhoSave[pt] + (-4*rhoSave[pt] + rhoSave[pt + unit0] + rhoSave[pt - unit0] +
                        + rhoSave[pt + unit1] + rhoSave[pt - unit1])/12);
    }
  MDWrite("residual",LRho);
};
int main(int argc, char* argv[])
{
  int M; 
  sscanf(argv[1],"%d",&M);
  // sscanf(argv[2],"%d",&inputMode);
  hockneyTest(M);
  return 0;
};
