CFLAGS += -I.. -I$(HOME)/local/src/atmosphere -std=c99 -O3 -g -Wall -MD

-include atmosphere.d
atmosphere: Makefile atmosphere.c
	$(CC) $(CFLAGS) -DOMEGA=0 atmosphere.c -o atmosphere -lgsl -lgslcblas -lm
