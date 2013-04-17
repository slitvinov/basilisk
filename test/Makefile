VALGRIND = valgrind -q --tool=memcheck
# the 'sse' stuff is a trick to turn off 'extended precision', 
# see: http://www.network-theory.co.uk/docs/gccintro/gccintro_70.html
CFLAGS += -std=c99 -g -O2 -MD -Wall -msse2 -mfpmath=sse

fast: halo.vlog speed.slog interp.log interpu.log refine.vlog rk2.vlog poisson.log circle.log \
	laplacian.slog interpolate.vlog events.vlog transcritical.log parabola.log bump.log \
	dry.log
	cat *.d > Makefile.deps

slow: advection.slog rotation.slog

plots: halo.png speed.png interp.png interpu.png refine.png rk2.png poisson.png circle.png \
	laplacian.png interpolate.png transcritical.png parabola.png dry.png

-include Makefile.deps

poisson.log: poisson.clog

laplacian.slog: laplacian.clog

clean:
	rm -f *.d *.log *.vlog *.slog *.clog *.out *.png \
	*.exe *.cexe *.pgexe *.prof Makefile.deps

%.png: %.plot
	gnuplot -e "batch=1; set term pngcairo enhanced; set output '$@';" $<

%.log: %.exe
	-rm -f $*.png
	./$< 2> $@ > $*.out
	diff $@ $*.ref

%.slog: %.exe
	-rm -f $*.png
	./$< 2> $@ > $*.out
	-diff $@ $*.ref

%.clog: %.cexe
	./$< 2> $@ > $*.cout

%.vlog: %.exe
	-rm -f $*.png
	$(VALGRIND) ./$< 2> $@ > $*.out
	diff $@ $*.ref

%.exe: %.c ../qcc
	../qcc $(CFLAGS) -o $@ -lm $<

%.cexe: %.c ../qcc
	../qcc -grid=multigrid $(CFLAGS) -o $@ -lm $<

%.pgexe: %.c ../qcc
	../qcc $(CFLAGS) -pg -o $@ -lm $<

%.prof: %.pgexe
	./$< 2> $*.log > $*.out
	gprof $< > $@
