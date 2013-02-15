qcc: lex.yy.c
	$(CC) $(CFLAGS) -DLIBDIR=\"`pwd`\" -O2 lex.yy.c -o qcc

lex.yy.c: qcc.lex
	flex qcc.lex
