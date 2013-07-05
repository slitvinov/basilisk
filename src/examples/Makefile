# the default CFLAGS are set in $(BASILISK)/config
CFLAGS += -g -O2 -fopenmp

all: check

tsunami.exe: tsunami.c $(BASILISK)/qcc Makefile
	qcc $(CFLAGS) tsunami.c -o tsunami.exe $(BASILISK)/kdt/kdt.o -lm

tides.exe: tides.c $(BASILISK)/qcc Makefile
	qcc $(CFLAGS) tides.c -o tides.exe $(BASILISK)/kdt/kdt.o -lm

include $(BASILISK)/Makefile.defs
