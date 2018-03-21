function shortline(line, string)
{
    i = match(line, string)
    j = i + length(string) - 60
    if (j > 0)
	line = "..." substr(line, j);
    if (length(line) > 80)
	line = substr(line, 1, 80) "...";
    return line;
}
    
BEGIN {
}
{
    file = $0
    sub(/:.*$/, "", file)
    sub(/^\./, "", file)
    line = $0
    sub(/^[^:]*:/, "", line)
    files[file] = file;
    lines[file][n[file]] = line;
    n[file]++;
}
END {
    N = 0
    for (file in files)
	N += n[file];
    print "### " N " matches found for <span id=pattern>" string "</span>";
    print "<ol>"
    for (file in files) {
	short = file
	sub(/\.page$/, "", short)
	veryshort = short
	sub(/^\//, "", veryshort)
	# This requires javascript to work
	# print "1. <a href=\"" short "\">" veryshort \
	#     "</a> (" n[file] " matching lines) " \
	#     "<a class=\"showmatch\" href=\"#\" style=\"\" "\
	#     "onclick=\"toggleMatches($(this));\">[show matches]</a>";
	# print "<pre class=matches style=\"display: none;\">"
	print "<li><a href=\"" short "\">" veryshort \
	    "</a> (" n[file] " matching lines)"
	print "<pre class=matches style=\"\">"
	for (line in lines[file]) {
	    highlighted = shortline(lines[file][line], string);
	    gsub(string, "<span class=highlighted>" string "</span>",
		 highlighted)
	    print highlighted
	}
	print "</pre></li>"
    }
    print "</ol>"
}
