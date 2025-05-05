# Vortex_Method

Implementation of the Vortex Method for simulating incompressible fluid flows.


The fftTools tests are almost the same as our homework 3 FFT implementation.

`make hockney` will build the Hockney version of the code. Then `./Hockney/testHockney2D.exe 8` will run the Hockney tests.

Use `testVortex_{x}.sh` to run the vortex method tests.
`x` (01 to 10) is the number of the test.

`make clean` or `make realclean` (more thorough) will remove all the object files and executables.



```
(base) Guos-MacBook-Pro-2:Vortex_Method guosiqi$ git diff --stat-width=80 --stat-graph-width=40  \
> --stat f982 1ef2 -- '*.h' '*.cpp' '*.wsd'
 Hockney/CutoffKernel.cpp            |  50 +++
 Hockney/Hockney.cpp                 | 206 ++++++------
 Hockney/HockneyTest.cpp             | 628 +++++++++++++++++++++++++++++++-----
 RectArray/DBox.cpp                  | 217 +++++++++++++
 Writers/WriteRectMDArray.cpp        | 147 ++++-----
 design/v0/Particle_in_Cell.wsd      |  66 ++++
 design/v0/Vortex_Particle.wsd       | 194 +++++++++++
 design/v1/system.wsd                | 230 +++++++++++++
 fftTools/FFT1DBRI.cpp               |  37 ---
 fftTools/FFT1DRecursive.cpp         |  76 -----
 fftTools/FFT1DTest.cpp              | 125 +++++++
 fftTools/FFTCTBRI.cpp               |  64 ----
 fftTools/FFTMD.cpp                  | 314 ++++++++----------
 fftTools/FFTMDTest.cpp              | 102 ++++++
 fftTools/FFTW1D.cpp                 |  53 ++-
 fftTools/PowerItoI.cpp              |  22 +-
 fftTools/tmp/PowerItoI.cpp          |  14 -
 vortexMethod/ParticleSet.cpp        |  37 +++
 vortexMethod/ParticleVelocities.cpp | 265 +++++++++++++++
 vortexMethod/VortexTest.cpp         | 548 +++++++++++++++++++++++--------
 20 files changed, 2599 insertions(+), 796 deletions(-)
 ```