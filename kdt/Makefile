include ../config

all: xyz2kdt kdtquery

kdt.o: kdt.c
	$(CC) $(CFLAGS) -D_FILE_OFFSET_BITS=64 -c kdt.c

xyz2kdt: xyz2kdt.c kdt.o kdt.h
	$(CC) $(CFLAGS) xyz2kdt.c kdt.o -o xyz2kdt -lm

kdtquery: kdtquery.c kdt.o kdt.h
	$(CC) $(CFLAGS) xyz2kdt.c kdt.o -o kdtquery -lm

clean:
	rm -f *.o
