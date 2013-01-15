CFLAGS += -I$(HOME)/local/src/atmosphere -std=c99 -O3 -g -Wall 

atmosphere: Makefile advection.c parameters.h init.h output.h grid.h
	$(CC) $(CFLAGS) -DOMEGA=0 advection.c -o advection -lm
