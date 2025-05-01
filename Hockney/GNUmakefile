# Makefile for Hockney directory

include ../Make.include

VPATH= . $(WRITER_HOME) $(RMDA_HOME) $(TIMER_HOME) $(FFT_HOME)

CFLAGS 		+= -w
CFLAGS 		+= -I. -I$(WRITER_HOME) -I$(TIMER_HOME) -I$(RMDA_HOME) -I$(FFT_HOME) -I$(VORTEX_HOME) -DDIM=$(DIM)

SRCFILES	:= $(wildcard $(TIMER_HOME)/*.cpp $(RMDA_HOME)/*.cpp $(WRITER_HOME)/*.cpp Hockney.cpp CutoffKernel.cpp)
OBJS 		:=$(patsubst %.cpp, %.o, $(SRCFILES)) 

EOBJS		:= $(addsuffix $(DIM)D.o, HockneyTest)


libfft%D.a:$(wildcard $(FFT_HOME)/*.H $(FFT_HOME)/*.cpp)
	@echo "Compiling libfft$(DIM)D.a"
	cd $(FFT_HOME);mkdir -p $(LIBS_LOCAL);make libfft$(DIM)D.a

libfft1D.a:$(wildcard $(FFT_HOME)/*.H $(FFT_HOME)/*.cpp)
	@echo "Compiling libfft1D.a"
	cd $(FFT_HOME);rm *.o;make libfft1D.a


all: hockneyTest

hockneyTest: GNUmakefile $(OBJS) HockneyTest$(osuffix) libfft$(DIM)D.a libfft1D.a
	@echo "Compiling HockneyTest"
	$(CXX) $(CFLAGS) $(FFTW_INCFLAGS) $(OBJS) $(LIBS) HockneyTest$(osuffix) -o testHockney$(DIM)D.exe

clean:
	rm -f *.o *.exe *.d *.vtk

realclean:
	rm -f *.o *.exe *.d $(LIBS_LOCAL)/*.a 
	cd $(FFT_HOME);make clean
	cd $(TIMER_HOME);make clean


-include $(OBJS:.o=.d)
-include $(EOBJS:.o=.d)



