# the default CFLAGS are set in $(BASILISK)/config
CFLAGS += -g -O2 -fopenmp

all: check

tsunami.tst: tsunami.s
	@LIBS=$(BASILISK)/kdt/kdt.o $(BASILISK)/dotest tsunami.tst

tides.tst: tides.s
	@LIBS=$(BASILISK)/kdt/kdt.o $(BASILISK)/dotest tides.tst

include $(BASILISK)/Makefile.defs
