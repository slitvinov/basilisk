CFLAGS += -O2 $(C99FLAGS) -D_POSIX_SOURCE -D_BSD_SOURCE

all: qcc qplot

qcc: qcc.c include.o
	$(CC) $(CFLAGS) -DCC=\"$(CC)\" -DC99FLAGS=\"$(C99FLAGS)\" \
		qcc.c include.o -o qcc

include.o: include.c
	$(CC) $(CFLAGS) -DLIBDIR=\"`pwd`\" -c include.c

qplot.o: qplot.c
	$(CC) $(CFLAGS) -c qplot.c

include.c: include.lex
	flex -P inc -o include.c include.lex

qcc.c: qcc.lex
	flex -o qcc.c qcc.lex

tags:
	etags *.h grid/*.h

dist:
	darcs dist
