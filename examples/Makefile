VALGRIND = valgrind -q --tool=memcheck
# the default CFLAGS are set in ../config
CFLAGS += -fpe -fopenmp `pwd`/../kdt/kdt.o -lm

all: check

include ../Makefile.defs
