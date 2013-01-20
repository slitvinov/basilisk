VALGRIND = valgrind -q --tool=memcheck
CFLAGS += -I$(HOME)/local/src/atmosphere -std=c99 -g -O2 -Wall -MD
LIBS = -lm

all: halo.log speed.log interp.log refine.log

halo.log: halo halo.plot
	$(VALGRIND) ./halo 2> halo.log > halo.out
	gnuplot -e "set term pngcairo; set output 'halo.png';" halo.plot
	diff halo.log halo.ref

-include halo.d
halo: halo.c
	$(CC) $(CFLAGS) halo.c -o halo -lm

speed.log: speed speed.plot
	./speed 2> speed.log > speed.out
	gnuplot -e "set term pngcairo; set output 'speed.png';" speed.plot
	-diff speed.log speed.ref

-include speed.d
speed: speed.c
	$(CC) $(CFLAGS) speed.c -o speed -lm

interp.log: interp interp.plot
	./interp 2> interp.log > interp.out
	gnuplot -e "set term pngcairo; set output 'interp.png';" interp.plot
	diff interp.log interp.ref

-include interp.d
interp: interp.c
	$(CC) $(CFLAGS) interp.c -o interp -lm

refine.log: refine refine.plot
	$(VALGRIND) ./refine 2> refine.log > refine.out
	gnuplot -e "set term pngcairo; set output 'refine.png';" refine.plot
	diff refine.log refine.ref

-include refine.d
refine: refine.c
	$(CC) $(CFLAGS) refine.c -o refine -lm

rk2.log: rk2 rk2.plot
	./rk2 2> rk2.log > rk2.out
	gnuplot -e "set term pngcairo; set output 'rk2.png';" rk2.plot
	diff rk2.log rk2.ref

-include rk2.d
rk2: rk2.c
	$(CC) $(CFLAGS) rk2.c -o rk2 -lm

test.log: test
	$(VALGRIND) ./test 2> test.log > test.out
#	gnuplot -e "set term pngcairo; set output 'test.png';" test.plot
#	diff test.log test.ref

-include test.d
test: test.c
	$(CC) $(CFLAGS) test.c -o test -lm



clean:
	rm -f *.d *.log *.out *.png \
	halo speed interp refine

