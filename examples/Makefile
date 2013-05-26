VALGRIND = valgrind -q --tool=memcheck
# the default CFLAGS are set in ../config
CFLAGS += -fpe

all: check

tsunami.exe: tsunami.c ../qcc Makefile
	../qcc $(CFLAGS) -o tsunami.exe tsunami.c `pwd`/../kdt/kdt.o -lm
tsunami.pgexe: tsunami.c ../qcc Makefile
	../qcc $(CFLAGS) -pg -o tsunami.pgexe tsunami.c `pwd`/../kdt/kdt.o -lm

include ../Makefile.defs
