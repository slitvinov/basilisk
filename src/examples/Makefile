# the default CFLAGS are set in $(BASILISK)/config
CFLAGS += -O2 -fopenmp

all: check

tsunami.tst: tsunami.s
	@LIBS="-L$(BASILISK)/kdt -lkdt" $(BASILISK)/runtest tsunami.tst

madsen-sv.c: madsen.c
	ln -s madsen.c madsen-sv.c 

madsen-sv.tst: madsen-sv.s
	@CFLAGS="-O2 -DSAINT_VENANT=1" $(BASILISK)/runtest madsen-sv.tst

madsen.tst: madsen.s madsen-sv.tst
	@CFLAGS="-O2" $(BASILISK)/runtest madsen.tst

include $(BASILISK)/Makefile.defs
