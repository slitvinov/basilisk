VALGRIND = valgrind -q --tool=memcheck
CFLAGS += -fpe -g -MD

# list only non-default (!.tst) tests
check: halo.vtst speed.stst refine.vtst rk2.vtst \
       laplacian.stst interpolate.vtst events.vtst \
       faces.vtst

slow: advection.tst rotation.tst bump2D.tst lid.tst

poisson.tst: poisson.ctst

laplacian.stst: laplacian.ctst

include Makefile.defs
