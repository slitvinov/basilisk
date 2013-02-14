CFLAGS += -std=c99 -O2 -g -Wall

atmosphere.log: atmosphere
	./atmosphere 2> atmosphere.log
	gnuplot -e "set term pngcairo; set output 'atmosphere.png';" atmosphere.plot
	-diff atmosphere.log atmosphere.ref

-include atmosphere.d
atmosphere: atmosphere.c parameters.h
	../endfor -grid=cartesian $(CFLAGS) -DOMEGA=0. -o atmosphere atmosphere.c -lgsl -lgslcblas -lm
