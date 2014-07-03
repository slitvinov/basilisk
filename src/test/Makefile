# the default CFLAGS are set in $(BASILISK)/config
CFLAGS += -O2

all: check slow

# list only non-default (!.tst) tests
check: halo.vtst rk2.vtst interpolate.vtst events.vtst faces.vtst

# these tests will not be done by the default 'check' target
slow: nonlinear.tst advection.tst rotation.tst \
	bump2D.tst hump.tst drybump2D.tst ponds.tst conical.tst

# these tests have special depencies/compilation requirements

seawallsv.c: seawall.c
	ln -s seawall.c seawallsv.c 

seawallsv.tst: seawallsv.s
	@CFLAGS="$(CFLAGS) -DSAINT_VENANT=1" $(BASILISK)/runtest seawallsv.tst

seawall.tst: seawallsv.tst

bore.tst: bore1.tst

explosion.tst: explosion.ctst

laplacian.tst: laplacian.ctst

lidmac.c: lid.c
	ln -s lid.c lidmac.c

lidmac.tst: lidmac.s
	@CFLAGS="$(CFLAGS) -DMAC=1" $(BASILISK)/runtest lidmac.tst

lid.tst: lidmac.tst

nonlinear.tst: nonlinear.s
	@LIBS="-lgsl -lgslcblas" $(BASILISK)/runtest nonlinear.tst

ponds.tst: ponds.s
	@LIBS="-L$(BASILISK)/kdt -lkdt" $(BASILISK)/runtest ponds.tst

poisson.tst: poisson.ctst

reversed.tst: reversed.ctst

rotate.tst: rotate.ctst

terrain.tst: terrain.s
	@LIBS="-L$(BASILISK)/kdt -lkdt" $(BASILISK)/runtest terrain.tst

include $(BASILISK)/Makefile.defs
