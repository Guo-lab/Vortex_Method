# Vortex_Method

Implementation of the Vortex Method for simulating incompressible fluid flows.


```
(base) Guos-MacBook-Pro-2:Vortex_Method guosiqi$ git diff --stat-width=80 --stat-graph-width=40  --stat f982 2c29 -- '*.h' '*.cpp'
 Hockney/CutoffKernel.cpp            |  50 +++++
 Hockney/Hockney.cpp                 | 206 ++++++++++---------
 Hockney/HockneyTest.cpp             | 398 +++++++++++++++++++++++++++++-------
 RectArray/DBox.cpp                  | 217 ++++++++++++++++++++
 Writers/WriteRectMDArray.cpp        | 147 ++++++-------
 fftTools/FFT1DBRI.cpp               |  37 ----
 fftTools/FFT1DRecursive.cpp         |  76 -------
 fftTools/FFT1DTest.cpp              | 125 +++++++++++
 fftTools/FFTCTBRI.cpp               |  64 ------
 fftTools/FFTMD.cpp                  | 314 +++++++++++++---------------
 fftTools/FFTMDTest.cpp              | 102 +++++++++
 fftTools/FFTW1D.cpp                 |  53 ++---
 fftTools/PowerItoI.cpp              |  22 +-
 fftTools/tmp/PowerItoI.cpp          |  14 --
 vortexMethod/ParticleSet.cpp        |  37 ++++
 vortexMethod/ParticleVelocities.cpp | 265 ++++++++++++++++++++++++
 vortexMethod/VortexTest.cpp         | 397 ++++++++++++++++++++++-------------
 17 files changed, 1728 insertions(+), 796 deletions(-)
 ```