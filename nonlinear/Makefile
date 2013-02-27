CFLAGS += -std=c99 -O2 -g -Wall -MD

nonlinear.log: nonlinear
	./nonlinear 2> nonlinear.log
	-diff nonlinear.log nonlinear.ref
	gnuplot -e "set term pngcairo; set output 'nonlinear.png';" nonlinear.plot

-include nonlinear.d
nonlinear: nonlinear.c ../qcc
	../qcc -grid=cartesian $(CFLAGS) -DOMEGA=0. -o nonlinear nonlinear.c -lgsl -lgslcblas -lm

clean:
	rm -f nonlinear
