# the default CFLAGS are set in ../config
CFLAGS += -fpe

all: tsunami.tst

tsunami.exe: tsunami.c ../qcc Makefile
	../qcc $(CFLAGS) -o tsunami.exe tsunami.c `pwd`/../kdt/kdt.o -lm

tsunami.tst: tsunami.exe
	./tsunami.exe > tsunami.out 2> tsunami.log

tsunami.png: tsunami.plot
	gnuplot -e "batch=1; set term pngcairo enhanced; set output '$@';" $< || rm -f $@
