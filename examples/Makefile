VALGRIND = valgrind -q --tool=memcheck
# the default CFLAGS are set in ../config
CFLAGS += -fpe -fopenmp

all: check

tsunami.exe: tsunami.c ../qcc Makefile
	../qcc $(CFLAGS) tsunami.c -o tsunami.exe `pwd`/../kdt/kdt.o -lm

tides.exe: tides.c ../qcc Makefile
	../qcc $(CFLAGS) tides.c -o tides.exe `pwd`/../kdt/kdt.o -lm

include ../Makefile.defs
