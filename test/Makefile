VALGRIND = valgrind -q --tool=memcheck
# the 'sse' stuff is a trick to turn off 'extended precision', 
# see: http://www.network-theory.co.uk/docs/gccintro/gccintro_70.html
CFLAGS += -std=c99 -g -O2 -Wall -msse2 -mfpmath=sse

all: halo.vlog speed.slog interp.log refine.vlog rk2.vlog poisson.log circle.log \
	laplacian.slog interpolate.vlog

plots: halo.png speed.png interp.png refine.png rk2.png poisson.png circle.png \
	laplacian.png interpolate.png

-include halo.d
-include speed.d
-include interp.d
-include refine.d
-include rk2.d
-include poisson.d
-include multigrid.d
-include circle.d
-include laplacian.d
-include interpolate.d

poisson.log: poisson.clog

laplacian.slog: laplacian.clog

clean:
	rm -f *.d *.log *.vlog *.slog *.clog *.out *.png *.exe *.cexe *.pgexe *.prof

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

%.exe: %.c
	../endfor $(CFLAGS) -o $@ -lm $<

%.cexe: %.c
	../endfor -grid=multigrid $(CFLAGS) -o $@ -lm $<

%.pgexe: %.c
	../endfor $(CFLAGS) -pg -o $@ -lm $<

%.prof: %.pgexe
	./$< 2> $*.log > $*.out
	gprof $< > $@
