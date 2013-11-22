# the default CFLAGS are set in $(BASILISK)/config
CFLAGS += -g -O2 -DTRASH=1

all: check slow

# list only non-default (!.tst) tests
check: halo.vtst rk2.vtst interpolate.vtst events.vtst faces.vtst

# these tests will not be done by the default 'check' target
slow: nonlinear.tst advection.tst rotation.tst \
	bump2D.tst hump.tst drybump2D.tst ponds.tst conical.tst

# these tests have special depencies/compilation requirements

laplacian.tst: laplacian.ctst

nonlinear.tst: nonlinear.s
	@LIBS="-lgsl -lgslcblas" $(BASILISK)/dotest nonlinear.tst

ponds.tst: ponds.s
	@LIBS=$(BASILISK)/kdt/kdt.o $(BASILISK)/dotest ponds.tst

poisson.tst: poisson.ctst

reversed.tst: reversed.ctst

rotate.tst: rotate.ctst

terrain.tst: terrain.s
	@LIBS=$(BASILISK)/kdt/kdt.o $(BASILISK)/dotest terrain.tst

include $(BASILISK)/Makefile.defs
