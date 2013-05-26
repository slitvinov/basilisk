VALGRIND = valgrind -q --tool=memcheck
# the default CFLAGS are set in ../config
CFLAGS += -fpe -DTRASH=1

all: check slow

# list only non-default (!.tst) tests
check: halo.vtst speed.stst refine.vtst rk2.vtst \
       laplacian.stst interpolate.vtst events.vtst \
       faces.vtst

slow: nonlinear.tst advection.tst rotation.tst \
	bump2D.tst hump.tst drybump2D.tst ponds.tst conical.tst

poisson.tst: poisson.ctst
laplacian.stst: laplacian.ctst
nonlinear.exe: nonlinear.c ../qcc Makefile
	../qcc $(CFLAGS) -o nonlinear.exe nonlinear.c -lgsl -lgslcblas -lm
terrain.exe: terrain.c ../qcc Makefile
	../qcc $(CFLAGS) -o terrain.exe terrain.c `pwd`/../kdt/kdt.o -lm
ponds.exe: ponds.c ../qcc Makefile
	../qcc $(CFLAGS) -o ponds.exe ponds.c `pwd`/../kdt/kdt.o -lm

include ../Makefile.defs
