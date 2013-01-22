all: endfor adapt.h atmosphere.h utils.h wavelet.h

.DELETE_ON_ERROR:

%.h: %.c endfor
	./endfor $< > $@
