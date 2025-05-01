# Makefile for vortexMethod directory

include ../Make.include

VPATH= . $(WRITER_HOME) $(RMDA_HOME) $(TIMER_HOME) $(FFT_HOME) $(HOCKNEY_HOME)

CFLAGS 		+= -I. -I$(WRITER_HOME) -I$(TIMER_HOME) -I$(RMDA_HOME) -I$(FFT_HOME) -I$(HOCKNEY_HOME) -DDIM=$(DIM)
CFLAGS 		+= -w

SRCFILES	:= $(wildcard $(TIMER_HOME)/*.cpp $(RMDA_HOME)/*.cpp $(WRITER_HOME)/*.cpp ./*.cpp)
SRCFILES 	+= $(HOCKNEY_HOME)/Hockney.cpp $(HOCKNEY_HOME)/CutoffKernel.cpp
OBJS		:=$(patsubst %.cpp,%.o, $(SRCFILES)) 


libfft$(DIM)D.a:$(wildcard $(FFT_HOME)/*.H $(FFT_HOME)/*.cpp)
	cd $(FFT_HOME);rm *.o;make libfft$(DIM)D.a

libfft1D.a:$(wildcard $(FFT_HOME)/*.H $(FFT_HOME)/*.cpp)
	cd $(FFT_HOME);rm *.o;make libfft1D.a

all: vortex2D

vortex2D: GNUmakefile libfft$(DIM)D.a libfft1D.a $(OBJS) 
	$(CXX) $(CFLAGS) $(FFTW_INCFLAGS) $(OBJS) $(LIBS) -o vortex$(DIM)D.exe

clean:
	rm -f *.o *.exe *.d

realclean:
	rm -f *.exe *.d *.o $(OBJS) $(LIBS_LOCAL)/*.a
	cd $(FFT_HOME);make clean
	cd $(TIMER_HOME);make clean
	cd $(HOCKNEY_HOME);make clean

-include $(OBJS:.o=.d)

