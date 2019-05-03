TEMPLATES=$(BASILISK)/darcsit/templates/

page.static: *.st $(TEMPLATES)includes.awk
	awk -v path=$(TEMPLATES) -f $(TEMPLATES)includes.awk page.st | \
	awk -v path=$(TEMPLATES) -f $(TEMPLATES)includes.awk | \
	awk -v path=$(TEMPLATES) -f $(TEMPLATES)includes.awk | \
	sed 's/\\\$$/__ESCAPEDDOLLAR__/g' > page.static
