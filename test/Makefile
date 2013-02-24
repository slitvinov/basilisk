VALGRIND = valgrind -q --tool=memcheck
# the 'sse' stuff is a trick to turn off 'extended precision', 
# see: http://www.network-theory.co.uk/docs/gccintro/gccintro_70.html
CFLAGS += -std=c99 -g -O2 -MD -Wall -msse2 -mfpmath=sse

all: halo.vlog speed.slog interp.log interpu.log refine.vlog rk2.vlog poisson.log circle.log \
	laplacian.slog interpolate.vlog events.vlog
	cat *.d > Makefile.deps

plots: halo.png speed.png interp.png interpu.png refine.png rk2.png poisson.png circle.png \
	laplacian.png interpolate.png

-include Makefile.deps

poisson.log: poisson.clog

laplacian.slog: laplacian.clog

clean:
	rm -f *.d *.log *.vlog *.slog *.clog *.out *.png \
	*.exe *.cexe *.pgexe *.prof Makefile.deps

%.png: %.plot
	gnuplot -e "batch=1; set term pngcairo; set output '$@';" $<

%.log: %.exe
	./$< 2> $@ > $*.out
	diff $@ $*.ref
	-rm -f $*.png

%.slog: %.exe
	./$< 2> $@ > $*.out
	-diff $@ $*.ref
	-rm -f $*.png

%.clog: %.cexe
	./$< 2> $@ > $*.cout

%.vlog: %.exe
	$(VALGRIND) ./$< 2> $@ > $*.out
	diff $@ $*.ref
	-rm -f $*.png

%.exe: %.c ../qcc
	../qcc $(CFLAGS) -o $@ -lm $<

%.cexe: %.c ../qcc
	../qcc -grid=multigrid $(CFLAGS) -o $@ -lm $<

%.pgexe: %.c ../qcc
	../qcc $(CFLAGS) -pg -o $@ -lm $<

%.prof: %.pgexe
	./$< 2> $*.log > $*.out
	gprof $< > $@
