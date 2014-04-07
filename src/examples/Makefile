# the default CFLAGS are set in $(BASILISK)/config
CFLAGS += -g -O2 -fopenmp

all: check

tsunami.tst: tsunami.c
	@LIBS="-L$(BASILISK)/kdt -lkdt" $(BASILISK)/runtest tsunami.tst

include $(BASILISK)/Makefile.defs
