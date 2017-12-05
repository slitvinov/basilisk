#/bin/sh

for f in `find . -name '*.[hc].page'` `find . -name '*.py.page'`; do
    l=`echo $f | sed 's/.page$//'`
    test -f $l || ln -f -s `basename $f` $l
done
