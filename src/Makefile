CFLAGS += -O2

# these are not tests
EXCLUDE = qcc.c include.c rotate.c qplot.c

all: links qcc qplot libkdt
	@chmod +x ppm2mpeg ppm2mp4 ppm2gif runtest sequence bview.page
	@test -f xyz2kdt || ln -s kdt/xyz2kdt
	@test -f kdtquery || ln -s kdt/kdtquery

bview_server: libgl bview-server2D bview-server3D bview2D bview3D

libkdt:
	cd kdt && make

libgl:
	cd gl && make

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

draw_get.h: draw.h params.awk
	awk -f params.awk < draw.h > draw_get.h

bview-server2D.d: bview-server.c qcc
	qcc -MD -o $@ $<

bview-server2D: bview-server.c bview-server2D.d draw_get.h
	qcc $(CFLAGS) bview-server.c -o $@ -Lgl -lglutils $(OPENGLIBS) -lm

bview-server3D.d: bview-server.c qcc
	qcc -grid=octree -MD -o $@ $<

bview-server3D: bview-server.c bview-server3D.d draw_get.h
	qcc -grid=octree $(CFLAGS) bview-server.c -o $@ \
		-Lgl -lglutils $(OPENGLIBS) -lm

bview2D:
	ln -s bview.page bview2D

bview3D:
	ln -s bview.page bview3D

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

-include bview-server3D.d
-include bview-server2D.d

include $(BASILISK)/Makefile.defs
