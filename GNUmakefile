# Makefile for the whole project

include Make.include

VPATH= . $(FFT_HOME) $(RDMA_HOME) $(LIBS_LOCAL) $(FFTW_HOME) $(TIMER_HOME) $(HOCKNEY_HOME) $(VORTEX_HOME) $(VTK_HOME)


.PHONY: timer
timer:
	cd $(TIMER_HOME);make clean;make all

.PHONY: fftTools
fftTools:
	cd $(FFT_HOME);make clean;make all

.PHONY: hockney
hockney:
	cd $(HOCKNEY_HOME);make clean;make all

.PHONY: vortex
vortex:
	cd $(VORTEX_HOME);make clean;make all



###############################
# Clean up
###############################
clean:
	rm -f *.o *.exe *.d *.vtk
	cd $(VTK_HOME); rm -f *.vtk

realclean:
	rm -f *.o *.exe *.d *.vtk
	cd $(TIMER_HOME);make clean
	cd $(FFT_HOME);make clean
	cd $(HOCKNEY_HOME);make clean
	cd $(VORTEX_HOME);make clean
	cd $(LIBS_LOCAL); rm -f *.a
	cd $(WRITER_HOME); rm -f *.o *.d
	cd $(RMDA_HOME); rm -f *.o *.d
	cd $(VTK_HOME); rm -f *.vtk
