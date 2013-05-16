# edit config (not this file) to tune compiler options etc..
include config

all: qcc qplot

qcc: qcc.c include.o
	$(CC99) $(CFLAGS) -DCC99="\"$(CC99)\"" qcc.c include.o -o qcc

include.o: include.c
	$(CC99) $(CFLAGS) -DLIBDIR=\"`pwd`\" -c include.c

qplot.o: qplot.c
	$(CC99) $(CFLAGS) -c qplot.c

include.c: include.lex
	flex -P inc -o include.c include.lex

qcc.c: qcc.lex
	flex -o qcc.c qcc.lex

tags:
	etags *.h grid/*.h

dist:
	darcs dist

clean:
	rm -f *.o
