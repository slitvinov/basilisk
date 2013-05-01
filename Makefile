all: qcc qplot

qcc: qcc.c include.o
	$(CC) $(CFLAGS) -O2 -DCFLAGS="\"$(CFLAGS)\"" qcc.c include.o -o qcc

include.o: include.c
	$(CC) $(CFLAGS) -O2 -DLIBDIR=\"`pwd`\" -c include.c

qplot.o: qplot.c
	$(CC) $(CFLAGS) -O2 -c qplot.c

include.c: include.lex
	flex -P inc -o include.c include.lex

qcc.c: qcc.lex
	flex -o qcc.c qcc.lex

tags:
	etags *.h grid/*.h
