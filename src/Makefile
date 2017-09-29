CFLAGS += -O2

# these are not tests
EXCLUDE = qcc.c include.c rotate.c qplot.c

all: links qcc qplot libkdt
	@chmod +x ppm2mpeg ppm2mp4 ppm2gif runtest sequence
	@test -f xyz2kdt || ln -s kdt/xyz2kdt
	@test -f kdtquery || ln -s kdt/kdtquery

libkdt:
	cd kdt && make

qcc: qcc.c include.o rotate.o config
	$(CC) $(CFLAGS) -DLIBDIR=\"`pwd`\" \
		-DCC99="\"$(CC99)\"" \
		-DCPP99="\"$(CPP99)\"" \
		-DBASILISK="\"$(BASILISK)\"" \
		qcc.c include.o rotate.o -o qcc

rotate.c: rotate.lex
	flex -P rot -o rotate.c rotate.lex

rotate: rotate.c
	$(CC) $(CFLAGS) -DTEST rotate.c -o rotate

rotate.o: rotate.c
	$(CC) $(CFLAGS) -c rotate.c

include.o: include.c
	$(CC) $(CFLAGS) -DLIBDIR=\"`pwd`\" -c include.c

qplot.o: qplot.c
	$(CC) $(CFLAGS) -c qplot.c

qplot: qplot.o
	$(CC) $(CFLAGS) qplot.o -o qplot

include.c: include.lex
	flex -P inc -o include.c include.lex

qcc.c: qcc.lex
	flex -o qcc.c qcc.lex

alltags:
	cd navier-stokes && make tags
	cd ehd && make tags
	cd examples && make tags
	cd test && make tags
	make tags
	cd navier-stokes && make itags
	cd ehd && make itags
	cd examples && make itags
	cd test && make itags
	make itags

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
