CFLAGS += -std=c99 -O2 -MD -g -Wall

lid.log: lid
	./lid 2> lid.log | ../qplot lid.qplot
	gnuplot -e "set term pngcairo; set output 'lid.png';" lid.plot
	-diff lid.log lid.ref

-include lid.d
lid: lid.c Makefile ../qcc
	../qcc -grid=multigrid $(CFLAGS) -o lid -lm lid.c

clean:
	rm -f *.log *.d
