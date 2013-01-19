CFLAGS += -I.. -I$(HOME)/local/src/atmosphere -std=c99 -O2 -g -Wall -MD

atmosphere.log: atmosphere
	./atmosphere 2> atmosphere.log
	gnuplot -e "set term pngcairo; set output 'atmosphere.png';" atmosphere.plot
	-diff atmosphere.log atmosphere.ref

-include atmosphere.d
atmosphere: atmosphere.c
	$(CC) $(CFLAGS) -DOMEGA=0 atmosphere.c -o atmosphere -lgsl -lgslcblas -lm
