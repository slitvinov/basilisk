#!/bin/sh

if test ! -d ../../_darcs; then
    touch Makefile.tests
    echo "Makefile.tests not updated (not darcs-controlled)"
else
    echo "updating Makefile.tests"
    (
    echo "# Automatically generated using 'make Makefile.tests'"
    echo "# DO NOT EDIT, edit 'Makefile' instead"
    echo "ALLTESTS = \\"
    DIR=`basename $PWD`
    darcs show files | grep 'src/'$DIR'/.*\.c' | 		\
	sed 's/.\/src\/'$DIR'\/\(.*\)/\t\1 \\/g'
    echo ""
    echo "PLOTS = \\"
    darcs show files | grep 'src/'$DIR'/.*\.plot' | 		\
	sed 's/.\/src\/'$DIR'\/\(.*\)/\t\1 \\/g'
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
    ) > Makefile.tests
fi
