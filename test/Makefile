VALGRIND = valgrind -q --tool=memcheck
CFLAGS += -I$(HOME)/local/src/atmosphere -std=c99 -g -Wall -MD
LDFLAGS += -lm

halo.log: Makefile halo halo.plot
	$(VALGRIND) ./halo 2> halo.log
	gnuplot -e "set term pngcairo; set output 'halo.png';" halo.plot
	diff halo.log halo.ref

-include halo.d
halo: halo.c
	$(CC) $(CFLAGS) halo.c -o halo -lm

clean:
	rm -f *.d
