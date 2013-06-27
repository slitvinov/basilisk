# the default CFLAGS are set in $(BASILISK)/config
CFLAGS += -g -O2 -fpe -DTRASH=1

# note that you can separate compilation and running the tests by doing
# % make exe
# % make
all: check slow

# list only non-default (!.tst) tests
check: halo.vtst speed.stst rk2.vtst \
       laplacian.stst interpolate.vtst events.vtst \
       faces.vtst

# these tests will not be done by the default 'check' target
slow: nonlinear.tst advection.tst rotation.tst \
	bump2D.tst hump.tst drybump2D.tst ponds.tst conical.tst

# these tests have special depencies/compilation requirements
poisson.tst: poisson.ctst
laplacian.stst: laplacian.ctst
nonlinear.exe: nonlinear.c $(BASILISK)/qcc Makefile
	qcc $(CFLAGS) -o nonlinear.exe nonlinear.c -lgsl -lgslcblas -lm
terrain.exe: terrain.c $(BASILISK)/qcc Makefile
	qcc $(CFLAGS) -o terrain.exe terrain.c $(BASILISK)/kdt/kdt.o -lm
ponds.exe: ponds.c $(BASILISK)/qcc Makefile
	qcc $(CFLAGS) -o ponds.exe ponds.c $(BASILISK)/kdt/kdt.o -lm

include $(BASILISK)/Makefile.defs
