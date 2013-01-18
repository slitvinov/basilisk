CFLAGS += -I$(HOME)/local/src/atmosphere -std=c99 -O2 -g -Wall -MD

-include advection.d
advection: Makefile advection.c parameters.h init.h output.h
	$(CC) $(CFLAGS) advection.c -o advection -lm
