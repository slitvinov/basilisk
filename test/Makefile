VALGRIND = valgrind -q --tool=memcheck
CFLAGS += -I$(HOME)/local/src/atmosphere -std=c99 -g -O2 -Wall -MD

all: halo.vlog speed.slog interp.log refine.vlog rk2.vlog poisson.log

pictures: halo.png speed.png interp.png refine.png rk2.png poisson.png

-include halo.d
-include speed.d
-include interp.d
-include refine.d
-include rk2.d
-include poisson.d

poisson.log: poisson.exe multigrid.log
	./$< 9 2> $@ > $*.out
	diff $@ $*.ref
	-rm -f $*.png

multigrid.log: multigrid.exe
	./$< 9 2> $@ > $*.out
	diff $@ $*.ref
	-rm -f $*.png

poisson.exe: poisson.c
	$(CC) $(CFLAGS) -DGRID=\"quadtree.c\" $< -o $@ -lm

multigrid.exe: poisson.c
	$(CC) $(CFLAGS) -DGRID=\"multigrid.c\" $< -o $@ -lm

clean:
	rm -f *.d *.log *.vlog *.slog *.out *.png *.exe

.SUFFIXES: .log .vlog .slog .exe .png .plot

.plot.png:
	gnuplot -e "set term pngcairo; set output '$@';" $^

.exe.log:
	./$< 2> $@ > $*.out
	diff $@ $*.ref
	-rm -f $*.png

.exe.slog:
	./$< 2> $@ > $*.out
	-diff $@ $*.ref
	-rm -f $*.png

.exe.vlog:
	$(VALGRIND) ./$< 2> $@ > $*.out
	diff $@ $*.ref
	-rm -f $*.png

.c.exe:
	$(CC) $(CFLAGS) $< -o $@ -lm
