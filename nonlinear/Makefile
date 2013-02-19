CFLAGS += -std=c99 -O2 -g -Wall

nonlinear.log: nonlinear
	./nonlinear 2> nonlinear.log
	gnuplot -e "set term pngcairo; set output 'nonlinear.png';" nonlinear.plot
	-diff nonlinear.log nonlinear.ref

-include nonlinear.d
nonlinear: nonlinear.c ../qcc
	../qcc -grid=cartesian $(CFLAGS) -DOMEGA=0. -o nonlinear nonlinear.c -lgsl -lgslcblas -lm

clean:
	rm -f nonlinear
