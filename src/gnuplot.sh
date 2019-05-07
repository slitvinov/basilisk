SVG="svg enhanced font ',11'"
if ! gnuplot -e "batch=1; PNG=\"$PNG\"; SVG=\"$SVG\"; set macros; set term $SVG;" \
     plots > /dev/null 2> gnuplot.log; then
    test=$1
    if test -z "$test"; then
	test=`basename $PWD`
    fi
    gawk -v test="$test" '
    /line [0-9]+:/ {
      match ($0, "line ([0-9]+):(.*)", a);
      print test ".c:" a[1] ": warning: gnuplot:" a[2];
    }' < gnuplot.log
    rm -f gnuplot.log `find . -name '_plot*.svg' -size 0`
    exit 1;
fi
rm -f gnuplot.log `find . -name '_plot*.svg' -size 0`
