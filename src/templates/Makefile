# gitit data needs to be installed in /usr/share/gitit/data/templates/
# This path can be changed in includes.awk.

page.static: page.st includes.awk Makefile
	awk -f includes.awk page.st | \
	awk -f includes.awk | awk -f includes.awk | \
	sed 's/\\\$$/__ESCAPEDDOLLAR__/g' > page.static
