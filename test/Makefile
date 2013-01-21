VALGRIND = valgrind -q --tool=memcheck
CFLAGS += -I$(HOME)/local/src/atmosphere -std=c99 -g -O2 -Wall -MD

all: halo.vlog speed.slog interp.log refine.vlog rk2.vlog

pictures: halo.png speed.png interp.png refine.png rk2.png

-include halo.d
-include speed.d
-include interp.d
-include refine.d
-include rk2.d

clean:
	rm -f *.d *.log *.vlog *.slog *.out *.png *.exe

.SUFFIXES: .log .vlog .slog .exe .png

.log.png:
	gnuplot -e "set term pngcairo; set output '$@';" $*.plot
.slog.png:
	gnuplot -e "set term pngcairo; set output '$@';" $*.plot
.vlog.png:
	gnuplot -e "set term pngcairo; set output '$@';" $*.plot

.exe.log:
	./$< 2> $@ > $*.out
	diff $@ $*.ref

.exe.slog:
	./$< 2> $@ > $*.out
	-diff $@ $*.ref

.exe.vlog:
	$(VALGRIND) ./$< 2> $@ > $*.out
	diff $@ $*.ref

.c.exe:
	$(CC) $(CFLAGS) $< -o $@ -lm
