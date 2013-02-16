qcc: qcc.c include.o
	$(CC) $(CFLAGS) -O2 qcc.c include.o -o qcc

include.o: include.c
	$(CC) $(CFLAGS) -DLIBDIR=\"`pwd`\" -O2 -c include.c

include.c: include.lex
	flex -P inc -o include.c include.lex

qcc.c: qcc.lex
	flex -o qcc.c qcc.lex
