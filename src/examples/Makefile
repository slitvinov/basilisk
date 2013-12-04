# the default CFLAGS are set in $(BASILISK)/config
CFLAGS += -g -O2 -fopenmp

all: check

tsunami.tst: tsunami.s
	@LIBS=$(BASILISK)/kdt/kdt.o $(BASILISK)/runtest tsunami.tst

tides.tst: tides.s
	@LIBS=$(BASILISK)/kdt/kdt.o $(BASILISK)/runtest tides.tst

include $(BASILISK)/Makefile.defs
