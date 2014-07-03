# parser for gnuplot commands in Literate C pages

BEGIN { nplots = 0 }

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
    print "set output '_plot" nplots++ ".png'";
}
