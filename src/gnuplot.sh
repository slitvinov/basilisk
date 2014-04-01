# parser for gnuplot commands in Literate C pages

awk '
/^[ \t]*~~~/ {
    if (gnuplot)
	gnuplot = 0;
}

{
    if (gnuplot)
	print $0;
    else
	print "#", $0;
}

/^[ \t]*~~~gnuplot/ {
    gnuplot = 1;
}' | ( \
    gnuplot -e "batch=1; PNG=\"$PNG\";"				\
            -e " set term $PNG enhanced font \",10\"; set macros;" - && \
    touch plots )
