include ../config

# if you want to use a memory-buffered KDT, add -DBUFFERED_KDT=1 to CFLAGS
CFLAGS += -O2

all: xyz2kdt kdtquery

kdt.o: kdt.c
	$(CC99) $(CFLAGS) -D_FILE_OFFSET_BITS=64 -c kdt.c

xyz2kdt: xyz2kdt.c kdt.o kdt.h
	$(CC99) $(CFLAGS) xyz2kdt.c kdt.o -o xyz2kdt -lm

kdtquery: kdtquery.c kdt.o kdt.h
	$(CC99) $(CFLAGS) xyz2kdt.c kdt.o -o kdtquery -lm

clean:
	rm -f *.o
