CFLAGS += -std=c99 -O2 -g -Wall \
	-DGRID=\"grid/multigrid.h\" \
	-DSIM=\"lid.h\"

lid.log: lid
	./lid 2> lid.log
	gnuplot -e "set term pngcairo; set output 'lid.png';" lid.plot
	-diff lid.log lid.ref

-include lid.d
lid: lid.h
	../qcc $(CFLAGS) -o lid -lm navier-stokes.c
