# edit config (not this file) to tune compiler options etc..
include config

CFLAGS += -O2

all: qcc qplot libkdt

libkdt:
	cd kdt && make

qcc: qcc.c include.o
	$(CC99) $(CFLAGS) -DLIBDIR=\"`pwd`\" -DCC99="\"$(CC99)\"" \
		qcc.c include.o -o qcc

include.o: include.c
	$(CC99) $(CFLAGS) -DLIBDIR=\"`pwd`\" -c include.c

qplot.o: qplot.c
	$(CC99) $(CFLAGS) -c qplot.c

qplot: qplot.o
	$(CC99) $(CFLAGS) qplot.o -o qplot

include.c: include.lex
	flex -P inc -o include.c include.lex

qcc.c: qcc.lex
	flex -o qcc.c qcc.lex

tags:
	etags *.h grid/*.h

dist:
	darcs dist

diff:
	tar czvf diff.tgz `darcs whatsnew -s | awk '{print $$2}'`

clean:
	rm -f *.o
