# Makefile for fftTools directory

include ../Make.include

VPATH = $(TIMER_HOME) $(RMDA_HOME) $(WRITER_HOME)

C_FLAGS_TEST = $(CFLAGS)
C_FLAGS_TEST += -w -I. -I$(FFT_HOME) -I$(RMDA_HOME) -I$(TIMER_HOME) -I$(FFTW_HOME)/include

CFLAGS += -w -I. $(FFTW_INCFLAGS) -I$(TIMER_HOME) -I$(RMDA_HOME) -I$(WRITER_HOME)

all: libfft1D.a libfft$(DIM)D.a test1d test${DIM}d

ifneq ($(DIM), 1)
MDFFTOBJS = FFTMD$(osuffix)
libfft$(DIM)D.a: GNUmakefile $(MDFFTOBJS)
	$(LIBTOOL) libfft$(DIM)D.a $(MDFFTOBJS) 
	mkdir -p $(LIBS_LOCAL);mv libfft$(DIM)D.a $(LIBS_LOCAL)
endif

1DFFTOBJS = FFTW1D.o PowerItoI.o # FFT1DBRI.o FFTCTBRI.o FFT1DRecursive.o
libfft1D.a: GNUmakefile $(1DFFTOBJS)
	$(LIBTOOL) libfft1D.a $(1DFFTOBJS)
	mkdir -p $(LIBS_LOCAL);mv libfft1D.a $(LIBS_LOCAL)


LIB_FLAGS_TEST:= $(FFTW_LIBS) -L$(LIBS_LOCAL) -lfft1D
LIB_FLAGS_TEST_${DIM}D:= $(LIB_FLAGS_TEST) -lfft$(DIM)D

SRCFILES:= $(wildcard $(RMDA_HOME)/*.cpp $(TIMER_HOME)/*.cpp)
OBJS:=$(patsubst %.cpp, %.o, $(SRCFILES))

test1d: FFT1DTest.cpp GNUmakefile libfft1D.a
	$(CXX) $(C_FLAGS_TEST) FFT1DTest.cpp $(LIB_FLAGS_TEST) -o TEST1D.exe

test${DIM}d: GNUmakefile $(OBJS) FFTMDTest$(osuffix) libfft1D.a libfft$(DIM)D.a
	$(CXX) $(C_FLAGS_TEST) $(OBJS) FFTMDTest$(osuffix) $(LIB_FLAGS_TEST_${DIM}D) -o TEST${DIM}D.exe

clean:
	rm *.o *.d *.exe;rm -f $(LIBS_LOCAL)/*.a

# handle missing dependency files
-include $(MDFFTOBJS:.o=.d)
-include $(1DFFTOBJS:.o=.d)
