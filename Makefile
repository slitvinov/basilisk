atmosphere: atmosphere.c Makefile parameters.h init.h
	$(CC) -std=c99 -O3 -g -Wall atmosphere.c -o atmosphere -lm
