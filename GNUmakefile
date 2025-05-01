FFT_HOME = $(HOME)/fftTools
RMDA_HOME = $(HOME)/RectMDArray
LIBS_LOCAL = $(HOME)/lib
# You will need to set the FFTW_HOME variable for your system, include, and LIB_FLAGS as
# below to find fftw. 
FFTW_HOME = /Users/colella/fftw

VPATH= . $(FFT_HOME) $(RDMA_HOME) $(LIBS_LOCAL) $(FFTW_HOME)

system := $(shell uname)

BLASLIBFLAGS = -lcblas
ifeq ($(system),Darwin)
  BLASFLAGS = -framework Accelerate
endif
osuffix:=$(DIM)D.o
dsuffix:=$(DIM)D.d
#CFLAGS = -g
CFLAGS = -O3 
CFLAGS += -std=c++11 -I. -I$(FFT_HOME) -I$(RMDA_HOME) -I$(FFTW_HOME)/include -DDIM=$(DIM) 
ifeq ($(CXX) , clang++)
CFLAGS += -stdlib=libc++
endif

#LIB_FLAGS:= -L$(LIBS_LOCAL) -L$(FFTW_HOME)/lib -lfftw3 -lfft1D
LIB_FLAGS:= -L $(LIBS_LOCAL) -lfft1D

SRCFILES:= $(wildcard $(RMDA_HOME)/*.cpp)
OBJS:=$(patsubst %.cpp,%.o, $(SRCFILES))

EOBJS:= FFTTest1D.o $(addsuffix $(DIM)D.o, FFTMDTest)

%.o: %.cpp GNUmakefile
	$(CXX) -c $(CFLAGS) $< -o $@
	$(CXX) -MM $(CFLAGS) $< > $*.d

%$(osuffix): %.cpp GNUmakefile
	$(CXX) -c $(CFLAGS) $< -o $@
	$(CXX) -MM $(CFLAGS) $< > $*$(dsuffix)

ifneq ($(DIM),1)
$(LIBS_LOCAL)/libfft%D.a:$(wildcard $(FFT_HOME)/*.cpp)  $(wildcard $(FFT_HOME)/*.H) GNUmakefile $(HOME)/GNUmakefile $(FFT_HOME)/GNUmakefile $(RMDA_HOME)/*.H $(RMDA_HOME)/*.cpp
	cd $(FFT_HOME);make clean;make libfft$*D.a DIM=$* CXX=$(CXX)
endif

$(LIBS_LOCAL)/libfft1D.a:$(wildcard $(FFT_HOME)/*.cpp) $(wildcard $(FFT_HOME)/*.H) GNUmakefile $(HOME)/GNUmakefile $(FFT_HOME)/GNUmakefile
	cd $(FFT_HOME);make clean;make libfft1D DIM=1 CXX=$(CXX)
clean:
	rm  *.o *.exe *.d

realclean:
	rm  *.o *.exe *.d $(LIBS_LOCAL)/*.a $(FFT_HOME)/*.d $(RMDA_HOME)/*.o $(RMDA_HOME)/*.d
-include $(OBJS:.o=.d)
-include $(EOBJS:.o=.d)
