#include <cmath>
#include <complex>
#include <vector>
#include <cstdio>
#include <iostream>
using namespace std;
#include "PowerItoI.H"
#include "FFTCTBRI.H"

void FFTCTBRI(vector<complex<double> > & a_f,const int& a_M, const bool& a_isForward)
{
  int N = Power(2,a_M);
 
  vector<complex<double> > f(N);

  for (int j = 0; j < N;j++)
  {
    int jlevel = j;
    int factor = N/2;
    int indOfj;
    indOfj = 0;
    for (int level = 0; level < a_M; level++)
      {
        int i = jlevel%2;
        jlevel = jlevel/2;
        indOfj += factor*i;
        factor = factor/2;
      }
    f[indOfj] = a_f[j];
  }
  for (int j = 0; j < N; j+=2)
    {
      complex<double> fHat0 = f[j] + f[j+1];
      complex<double> fHat1 = f[j] - f[j+1];
      f[j] = fHat0;
      f[j+1] = fHat1;      
    }
  int levelLength = 2;
  int sign;
  if (a_isForward) {sign = -1;} else {sign=1;}
  for (int level = 1; level < a_M;level++)
    {
      double h = 1./levelLength/2;
      complex<double> z0Level(cos(2*M_PI*h),sign*sin(2*M_PI*h));
      for (int j0 = 0; j0 < N; j0+=2*levelLength)
        {
          complex<double> zToTheK(1.,0);
          for (int j = 0; j < levelLength; j++)
            {
              complex<double> temp = zToTheK*f[j + j0 + levelLength];
              complex<double> fHat0 = f[j + j0] + temp;
              complex<double> fHat1 = f[j + j0] - temp;
              f[j + j0] = fHat0;
              f[j + j0 + levelLength] = fHat1;
              zToTheK *= z0Level;
            }
        }
      levelLength *= 2;
    }
 for (int j = 0; j < N;j++)
  {
    a_f[j] = f[j];
  }
};
