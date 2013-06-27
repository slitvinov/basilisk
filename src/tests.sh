#!/bin/sh

showfiles()
{
    if ! darcs show files > /dev/null 2>&1; then
	ls *.$1 2> /dev/null | sed 's/\(.*\)/\t\1 \\/g'
    else
	DIR=`basename $PWD`
	darcs show files | grep '/'$DIR'/[^/]*\.'$1'$' | \
	    sed 's/.*\/'$DIR'\/\(.*\)/\t\1 \\/g'
    fi
}

echo "updating Makefile.tests"
(
    echo "# Automatically generated using 'make Makefile.tests'"
    echo "# DO NOT EDIT, edit 'Makefile' instead"
    DIR=`basename $PWD`
    echo "ALLTESTS = \\"
    showfiles c
    echo ""
    echo "PLOTS = \\"
    showfiles plot
    echo ""
    echo "TESTS = \\"
    sed ':x; /\\$/ { N; s/\\\n//; tx }' < Makefile | 		\
	grep '^check:' | grep -o '[a-zA-Z_0-9]*\..tst' | 	\
	sed 's/\(.*\)\..tst/\t\1.c \\/g'
    sed ':x; /\\$/ { N; s/\\\n//; tx }' < Makefile | 		\
	grep -v '^check:' | grep -o ':.*' | 			\
	grep -o '[a-zA-Z_0-9]*\.tst' | 				\
	sed 's/\(.*\)\.tst/\t\1.c \\/g'
    echo ""
    echo "SPECIAL_TESTS = \\"
    grep -o '[a-zA-Z_0-9]*\.ctst' Makefile | sort | uniq |	\
	sed 's/\(.*\)/\t\1 \\/g'
    echo ""
    echo "ALLPAGES = \\"
    showfiles '[ch]\.page'
    echo ""
) > Makefile.tests

rm -f Makefile.deps
