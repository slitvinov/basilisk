all: check 3D-tests mpi-tests curvature-tests axi-tests

include $(BASILISK)/Makefile.defs

# the default CFLAGS are set in $(BASILISK)/config
CFLAGS += -O2 -DMTRACE=3

# list only non-default (!.tst) tests
check: halo.vtst rk2.vtst interpolate.vtst events.vtst faces.vtst coarsen.vtst

# these tests have special depencies/compilation requirements

dry-explicit.c: dry.c
	ln -s dry.c dry-explicit.c 

dry-explicit.tst: dry-explicit.s
	@CFLAGS="$(CFLAGS) -DEXPLICIT=1" $(BASILISK)/runtest dry-explicit.tst

dry.tst: dry-explicit.tst

drybump-explicit.c: drybump.c
	ln -s drybump.c drybump-explicit.c 

drybump-explicit.tst: drybump-explicit.s
	@CFLAGS="$(CFLAGS) -DEXPLICIT=1" $(BASILISK)/runtest drybump-explicit.tst

drybump.tst: drybump-explicit.tst

drybump2D-implicit.c: drybump2D.c
	ln -s drybump2D.c drybump2D-implicit.c 

drybump2D-implicit.tst: drybump2D-implicit.s
	@CFLAGS="$(CFLAGS) -DIMPLICIT=1" \
		$(BASILISK)/runtest drybump2D-implicit.tst

drybump2D.tst: drybump2D-implicit.tst

parabola-explicit.c: parabola.c
	ln -s parabola.c parabola-explicit.c 

parabola-explicit.tst: parabola-explicit.s
	@CFLAGS="$(CFLAGS) -DEXPLICIT=1" $(BASILISK)/runtest \
		parabola-explicit.tst

parabola.tst: parabola-explicit.tst

seawallsv.c: seawall.c
	ln -s seawall.c seawallsv.c 

seawallsv.tst: seawallsv.s
	@CFLAGS="$(CFLAGS) -DSAINT_VENANT=1" $(BASILISK)/runtest seawallsv.tst

seawall.tst: seawallsv.tst

conicalsv.c: conical.c
	ln -s conical.c conicalsv.c 

conicalsv.tst: conicalsv.s
	@CFLAGS="$(CFLAGS) -DSAINT_VENANT=1" $(BASILISK)/runtest conicalsv.tst

conical-implicit.c: conical.c
	ln -s conical.c conical-implicit.c 

conical-implicit.tst: conical-implicit.s
	@CFLAGS="$(CFLAGS) -DSAINT_VENANT=1 -DIMPLICIT=1" \
		$(BASILISK)/runtest conical-implicit.tst

conical.tst: conicalsv.tst conical-implicit.tst

curvature.tst: curvature.s curvature.3D.tst
	@LIBS="-L$(HOME)/local/lib -lvofi" \
	$(BASILISK)/runtest curvature.tst

curvature.3D.tst: curvature.3D.s
	@CFLAGS="$(CFLAGS) -grid=octree" \
	LIBS="-L$(HOME)/local/lib -lvofi" \
	$(BASILISK)/runtest curvature.3D.tst

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

ponds-implicit.c: ponds.c
	ln -s ponds.c ponds-implicit.c

ponds-implicit.tst: ponds-implicit.s
	@LIBS="-L$(BASILISK)/kdt -lkdt" CFLAGS="$(CFLAGS) -DIMPLICIT=1" \
		$(BASILISK)/runtest ponds-implicit.tst

ponds.tst: ponds.s ponds-implicit.tst
	@LIBS="-L$(BASILISK)/kdt -lkdt" $(BASILISK)/runtest ponds.tst

poisson.tst: poisson.ctst

rising-axi.c: rising.c
	ln -s rising.c rising-axi.c

rising-axi.tst: rising-axi.s rising.s
	@CFLAGS="$(CFLAGS) -DAXIS=1" $(BASILISK)/runtest rising-axi.tst

reversed.tst: reversed.ctst

rotate.tst: rotate.ctst

taylor-green-all-mach.c: taylor-green.c
	ln -s taylor-green.c taylor-green-all-mach.c 

taylor-green-all-mach.tst: taylor-green-all-mach.s
	@CFLAGS="$(CFLAGS) -DALL_MACH=1" $(BASILISK)/runtest taylor-green-all-mach.tst

taylor-green.tst: taylor-green-all-mach.tst

terrain.tst: terrain.s
	@LIBS="-L$(BASILISK)/kdt -lkdt" $(BASILISK)/runtest terrain.tst

coarsen.tst: CFLAGS=-DMTRACE=3

# Axisymmetric tests

axi-tests: axiadvection.tst axi.tst plateau.tst poiseuille-axi.tst \
	rising-axi.tst

# 3D tests

3D-tests: circle.3D.tst curvature.3D.tst hf.3D.tst periodic.3D.tst \
	poisson.3D.tst refineu.3D.tst solenoidal.3D.tst

# Curvature tests

curvature-tests: hf.tst hf3D.tst hf.ctst hf.3D.tst hf1.tst \
	curvature.tst curvature.3D.tst

# MPI tests

mpi-tests: indexing.tst indexing.3D.tst \
	mpi-restriction.tst mpi-restriction.3D.tst \
	mpi-reduce.tst mpi-refine.tst mpi-refine.3D.tst \
	mpi-laplacian.tst mpi-laplacian.3D.tst \
	mpi-circle.tst mpi-circle1.tst mpi-flux.tst \
	mpi-interpu.tst mpi-coarsen.tst mpi-coarsen1.tst \
	bump2Dp.tst axiadvection.tst hf1.tst vortex.tst pdump.tst restore.tst \
	gfsi.tst gfs.tst

indexing.tst:		CC = mpicc -D_MPI=3
indexing.3D.tst:	CC = mpicc -D_MPI=3
mpi-restriction.tst:	CC = mpicc -D_MPI=3
mpi-restriction.3D.tst:	CC = mpicc -D_MPI=3
mpi-reduce.tst:		CC = mpicc -D_MPI=3
mpi-refine.tst:		CC = mpicc -D_MPI=4
mpi-refine.3D.tst:	CC = mpicc -D_MPI=4
mpi-laplacian.tst:	CC = mpicc -D_MPI=3
mpi-laplacian.3D.tst:	CC = mpicc -D_MPI=3
mpi-circle.tst:		CC = mpicc -D_MPI=5

mpi-circle1.c: mpi-circle.c
	ln -s mpi-circle.c mpi-circle1.c
bump2Dp.c: bump2D.c
	ln -s bump2D.c bump2Dp.c

mpi-circle1.tst:  CC = mpicc -D_MPI=6
mpi-flux.tst:     CC = mpicc -D_MPI=6
mpi-interpu.tst:  CC = mpicc -D_MPI=5
mpi-coarsen.tst:  CC = mpicc -D_MPI=2
mpi-coarsen1.tst: CC = mpicc -D_MPI=5
bump2Dp.tst:      CC = mpicc -D_MPI=9
vortex.tst:	  CC = mpicc -D_MPI=7
axiadvection.tst: CC = mpicc -D_MPI=7
hf1.tst:	  CC = mpicc -D_MPI=7

pdump.tst: CC = mpicc -D_MPI=7
restore.tst: pdump.tst

gfs.tst: CC = mpicc -D_MPI=7
gfsi.tst: gfs.tst
