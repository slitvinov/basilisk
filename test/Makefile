VALGRIND = valgrind -q --tool=memcheck
CFLAGS += -I$(HOME)/local/src/atmosphere -std=c99 -g -O2 -Wall -MMD

all: halo.vlog speed.slog interp.log refine.vlog rk2.vlog poisson.log

plots: halo.png speed.png interp.png refine.png rk2.png poisson.png

-include halo.d
-include speed.d
-include interp.d
-include refine.d
-include rk2.d
-include poisson.d

poisson.log: multigrid.log

poisson.exe: poisson.c
	../endfor poisson.c | $(CC) $(CFLAGS) -DGRID=\"grid/quadtree.h\" -x c - -o poisson.exe -lm

multigrid.exe: poisson.c
	../endfor poisson.c | $(CC) $(CFLAGS) -DGRID=\"grid/multigrid.h\" -x c - -o multigrid.exe -lm

clean:
	rm -f *.d *.log *.vlog *.slog *.out *.png *.exe

%.png: %.plot
	gnuplot -e "set term pngcairo; set output '$@';" $^

%.log: %.exe
	./$< 2> $@ > $*.out
	diff $@ $*.ref
	-rm -f $*.png

%.slog: %.exe
	./$< 2> $@ > $*.out
	-diff $@ $*.ref
	-rm -f $*.png

%.vlog: %.exe
	$(VALGRIND) ./$< 2> $@ > $*.out
	diff $@ $*.ref
	-rm -f $*.png

%.exe: %.c
	../endfor $< | $(CC) $(CFLAGS) -x c - -o $@ -lm
