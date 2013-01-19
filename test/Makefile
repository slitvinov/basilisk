VALGRIND = valgrind -q --tool=memcheck
CFLAGS += -I$(HOME)/local/src/atmosphere -std=c99 -g -O2 -Wall -MD
LIBS = -lm

all: halo.log speed.log interp.log

halo.log: halo halo.plot
	$(VALGRIND) ./halo 2> halo.log
	gnuplot -e "set term pngcairo; set output 'halo.png';" halo.plot
	diff halo.log halo.ref
-include halo.d
halo: halo.c
	$(CC) $(CFLAGS) halo.c -o halo -lm

speed.log: speed speed.plot
	./speed 2> speed.log
	gnuplot -e "set term pngcairo; set output 'speed.png';" speed.plot
	-diff speed.log speed.ref
-include speed.d
speed: speed.c
	$(CC) $(CFLAGS) speed.c -o speed -lm

interp.log: interp interp.plot
	$(VALGRIND) ./interp 2> interp.log
	gnuplot -e "set term pngcairo; set output 'interp.png';" interp.plot
	diff interp.log interp.ref
-include interp.d
interp: interp.c
	$(CC) $(CFLAGS) interp.c -o interp -lm

clean:
	rm -f *.d
