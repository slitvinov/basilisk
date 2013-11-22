CFLAGS += -O2

# these are not tests
EXCLUDE = qcc.c include.c rotate.c qplot.c

all: qcc qplot libkdt
	@chmod +x ppm2mpeg runtest
	@test -f xyz2kdt || ln -s kdt/xyz2kdt

libkdt:
	cd kdt && make

qcc: qcc.c include.o rotate.o
	$(CC99) $(CFLAGS) -DLIBDIR=\"`pwd`\" \
		-DCC99="\"$(CC99)\"" \
		-DCPP99="\"$(CPP99)\"" \
		qcc.c include.o rotate.o -o qcc

rotate.c: rotate.lex
	flex -P rot -o rotate.c rotate.lex

rotate: rotate.c
	$(CC99) $(CFLAGS) -DTEST rotate.c -o rotate

rotate.o: rotate.c
	$(CC99) $(CFLAGS) -c rotate.c

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

alltags: tags
	cd navier-stokes && make tags
	cd examples && make tags
	cd test && make tags

etags:
	etags *.h grid/*.h

changelog:
	darcs changes > ChangeLog

dist:
	darcs dist

diff:
	cd .. && tar czvf src/diff.tgz `darcs whatsnew -s | \
		sed 's/. .\/.*\/$$//g' | awk '{print $$2}'`

include $(BASILISK)/Makefile.defs
