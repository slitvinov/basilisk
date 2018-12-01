if ! gnuplot -e "batch=1; PNG=\"$PNG\"; SVG=\"svg enhanced font \\\",14\\\"\"; set term SVG; set macros;" \
    plots > /dev/null 2> gnuplot.log; then
    gawk -v test=`basename $PWD` '
    /line [0-9]+:/ {
      match ($0, "line ([0-9]+):(.*)", a);
      print test ".c:" a[1] ": warning: gnuplot:" a[2];
    }' < gnuplot.log
    rm -f gnuplot.log
    exit 1;
fi
rm -f gnuplot.log
