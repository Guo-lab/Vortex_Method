# Makefile for timer directory

include ../Make.include

CFLAGS += -w

help:
	@echo ""
	@echo "Makefile for timer directory"
	@echo "make CH_Timer.o to compile the CH_Timer object file"
	@echo "make clean to remove all object files and executables"
	@echo "make all to check if all object files are up to date"
	@echo ""
	@echo ""


all: CH_Timer.o
	@echo "All object files are up to date."

CH_Timer.o: CH_Timer.H CH_Timer.cpp ClockTicks.H GNUmakefile
	$(CXX) $(CFLAGS) -c -o CH_Timer.o CH_Timer.cpp

clean:
	rm -f *.o *.d
	@echo "Cleaned up object files and executables."