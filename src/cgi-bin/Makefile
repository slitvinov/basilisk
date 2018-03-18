# -*-Makefile-*-
# Only used for generating HTML listings (not pages)

%.html: $(BASILISK)/cgi-bin/Makefile 		\
	$(BASILISK)/templates/page.static 	\
	$(BASILISK)/page2html
	$(BASILISK)/page2html $*
