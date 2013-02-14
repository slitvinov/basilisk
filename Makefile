qcc: lex.yy.c
	$(CC) $(CFLAGS) -O2 lex.yy.c -o qcc

lex.yy.c: qcc.lex
	flex qcc.lex
