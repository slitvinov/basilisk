all: endfor

lex.yy.c: endfor.lex
	flex endfor.lex

endfor: lex.yy.c
	cc -O2 -Wall -Wno-unused -g lex.yy.c -o endfor -lfl
