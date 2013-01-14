quadtree: quadtree.c
	$(CC) -std=c99 -O3 -g -Wall quadtree.c -o quadtree -lm

atmosphere: atmosphere.c Makefile parameters.h init.h boundary.h
	$(CC) -std=c99 -O3 -g -Wall atmosphere.c -o atmosphere -lm
