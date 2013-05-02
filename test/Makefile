VALGRIND = valgrind -q --tool=memcheck
CFLAGS += -fpe -DTRASH=1 -g -Wall -O2

# list only non-default (!.tst) tests
check: halo.vtst speed.stst refine.vtst rk2.vtst \
       laplacian.stst interpolate.vtst events.vtst \
       faces.vtst

slow: nonlinear.tst advection.tst rotation.tst bump2D.tst

poisson.tst: poisson.ctst
laplacian.stst: laplacian.ctst
nonlinear.exe: nonlinear.c ../qcc
	../qcc $(CFLAGS) -o nonlinear.exe nonlinear.c -lgsl -lgslcblas -lm

include Makefile.defs
