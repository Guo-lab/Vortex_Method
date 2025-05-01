#our makefile variables.   Good to place commonly changed variables
# at the top of your makefile. Then follow with your rules.
VPATH= . ../Writers ../RectArray ../fftTools ../timer 
DIM=2
HOME = ..
system := $(shell uname)
FFTWDIR = /usr/local
FFTINCFLAGS = -I$(FFTWDIR)/include
CFLAGS = -g -Wall
osuffix:=$(DIM)D.o
dsuffix:=$(DIM)D.d
#CFLAGS = -O3 
CFLAGS += -I. -I../Writers -I../timer -I../RectArray -I../fftTools -I../vortexMethod -DDIM=$(DIM)  -std=c++11 
CFLAGS += $(FFTINCFLAGS)
LIBS:=  -L$(FFTWDIR)/lib -lfftw3 -L../lib -lfft$(DIM)D -lfft1D

#CXX = /Users/colella/Desktop/gcc-4.7-bin/usr/local/bin/g++ 
CXX=g++

SRCFILES:= $(wildcard ../timer/*.cpp ../RectArray/*.cpp ../Writers/*.cpp Hockney.cpp)
OBJS:=$(patsubst %.cpp,%.o, $(SRCFILES)) 

EOBJS:= $(addsuffix $(DIM)D.o, HockneyTest) 

%.o: %.cpp GNUmakefile
	$(CXX) -c $(CFLAGS) $< -o $@
	$(CXX) -MM $(CFLAGS) $< > $*.d

%$(osuffix): %.cpp GNUmakefile
	$(CXX) -c $(CFLAGS) $< -o $@
	$(CXX) -MM $(CFLAGS) $< > $*$(dsuffix)

../lib/libfft%D.a:$(wildcard ../fftTools/*.H ../fftTools/*.cpp)
	cd ../fftTools;mkdir -p ../lib;make DIM=$* CXX=$(CXX) FFTWDIR=$(FFTWDIR)

../lib/libfft1D.a:$(wildcard $(HOME)/fftTools/*.H $(HOME)/fftTools/*.cpp)
	cd $(HOME)/fftTools;rm *.o;make libfft1D.a CXX=$(CXX) DIM=$(DIM) FFTWDIR=$(FFTWDIR) FFTINCFLAGS=$(FFTINCFLAGS)

vars:
	echo $(wildcard ../fftTools/*.H ../fftTools/*.cpp)

hockneyTest: GNUmakefile $(OBJS) HockneyTest$(osuffix) ../lib/libfft$(DIM)D.a ../lib/libfft1D.a

	$(CXX) $(CFLAGS) $(FFTINCFLAGS) $(OBJS) $(LIBS) $(TOOLS) HockneyTest$(osuffix) -o testHockney$(DIM)D.exe

clean:
	rm  *.o *.exe *.d *.vtk

realclean:
	rm  *.o *.exe *.d ../lib/*.a ../fftTools/*.o ../fftTools/*.d ../lib/*

-include $(OBJS:.o=.d)
-include $(EOBJS:.o=.d)
