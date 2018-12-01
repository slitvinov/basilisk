# parser for gnuplot commands in Literate C pages

BEGIN { nplots = 0 }

function tightbb(output)
{
    if (0) # This is too expensive and requires inkscape
	print "! sed 's|^<rect x=\"0\" y=\"0\" width=\".*/>||' "	\
	    output							\
	    " | inkscape -z -D -l "					\
	    output							\
	    " -f /dev/stdin 2> /dev/null";
    else
	print "# ";
}

function defaults()
{
    printf "set pointsize 0.75; ";
    print "set key spacing 0.8";
}

/^[ \t]*~~~/ {
    if (gnuplot) {
	gnuplot = 0;
	if (match (output, ".*\\.png"))
	    print "! mogrify -trim" output;
	else if (match (output, ".*\\.svg"))
	    tightbb(output);
	else if (output == "")
	    tightbb("_plot" nplots++ ".svg");
    }
}

{
    if (gnuplot) {
	if (match($0, "^[ \t]*reset")) {
	    printf "reset; ";
	    printf "load '~/.gnuplot'; ";
	    defaults();
	}
	else {
	    if (match($0, "^[ \t]*set[ \t]+output[ \t]+(.*)", a))
		output = a[1];
	    print $0;
	}
    }
    else if (!match($0, "^[ \t]*~~~"))
	print "#", $0;
}

/^[ \t]*~~~gnuplot/ {
    gnuplot = 1; output = "";
    printf "set output '_plot" nplots ".svg'; ";
    defaults();
}
