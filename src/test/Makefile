all: check mpi-tests

include $(BASILISK)/Makefile.defs

# the default CFLAGS are set in $(BASILISK)/config
CFLAGS += -O2

# list only non-default (!.tst) tests
check: halo.vtst rk2.vtst interpolate.vtst events.vtst faces.vtst coarsen.vtst

# these tests have special depencies/compilation requirements

hf.tst: hf.ctst

seawallsv.c: seawall.c
	ln -s seawall.c seawallsv.c 

seawallsv.tst: seawallsv.s
	@CFLAGS="$(CFLAGS) -DSAINT_VENANT=1" $(BASILISK)/runtest seawallsv.tst

seawall.tst: seawallsv.tst

conicalsv.c: conical.c
	ln -s conical.c conicalsv.c 

conicalsv.tst: conicalsv.s
	@CFLAGS="$(CFLAGS) -DSAINT_VENANT=1" $(BASILISK)/runtest conicalsv.tst

conical.tst: conicalsv.tst

explicit.c: implicit.c
	ln -s implicit.c explicit.c 

explicit.tst: explicit.s
	@CFLAGS="$(CFLAGS) -DEXPLICIT=1" $(BASILISK)/runtest explicit.tst

implicit.tst: explicit.tst

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

rising-axi.c: rising.c
	ln -s rising.c rising-axi.c

rising-axi.tst: rising-axi.s
	@CFLAGS="$(CFLAGS) -DAXIS=1" $(BASILISK)/runtest rising-axi.tst

rising.tst: rising-axi.tst

reversed.tst: reversed.ctst

rotate.tst: rotate.ctst

terrain.tst: terrain.s
	@LIBS="-L$(BASILISK)/kdt -lkdt" $(BASILISK)/runtest terrain.tst

# MPI tests

mpi-tests: indexing.tst mpi-restriction.tst mpi-reduce.tst mpi-refine.tst \
	mpi-laplacian.tst mpi-circle.tst mpi-circle1.tst mpi-flux.tst \
	mpi-interpu.tst mpi-coarsen.tst bump2Dp.tst vortex.tst \
	axiadvection.tst

indexing.tst:		CC = mpicc -D_MPI=3
mpi-restriction.tst:	CC = mpicc -D_MPI=3
mpi-reduce.tst:		CC = mpicc -D_MPI=3
mpi-refine.tst:		CC = mpicc -D_MPI=4
mpi-laplacian.tst:	CC = mpicc -D_MPI=3
mpi-circle.tst:		CC = mpicc -D_MPI=5

mpi-circle1.c: mpi-circle.c
	ln -s mpi-circle.c mpi-circle1.c
bump2Dp.c: bump2D.c
	ln -s bump2D.c bump2Dp.c

mpi-circle1.tst: CC = mpicc -D_MPI=6
mpi-flux.tst:    CC = mpicc -D_MPI=6
mpi-interpu.tst: CC = mpicc -D_MPI=5
mpi-coarsen.tst: CC = mpicc -D_MPI=2
bump2Dp.tst:     CC = mpicc -D_MPI=9
vortex.tst:	 CC = mpicc -D_MPI=7
axiadvection.tst: CC = mpicc -D_MPI=7
