# parser for gnuplot commands in Literate C pages

BEGIN { nplots = 0 }

/^[ \t]*~~~/ {
    if (gnuplot)
	gnuplot = 0;
}

/^[ \t]*set[ \t]+output/ {
    if (gnuplot)
	output = 1;
}

/^[ \t]*(plot|splot)/ {
    if (gnuplot && !output)
	print "set output '_plot" nplots++ ".png'";
}

{
    if (gnuplot)
	print $0;
    else
	print "#", $0;
}

/^[ \t]*~~~gnuplot/ {
    gnuplot = 1;
    output = 0;
}
