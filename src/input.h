/**
# Importing Portable Gray Map (PGM) images

This function reads in the
[PGM](http://en.wikipedia.org/wiki/Netpbm_format) file in file *fp*
and imports the corresponding values into field *s*. The grayscale is
normalised and inverted so that the maximum value in field *s* is one
(black) and the minimum value is zero (white). 

By default the origin of the image (lower-left corner) is assumed to
be at (0,0) and the width of the image is set to `L0`. This can be
changed using the optional *ox*, *oy* and *width* parameters. */

struct InputPGM {
  // compulsory
  scalar s;
  FILE * fp;
  // optional
  double ox, oy, width;
};

void input_pgm (struct InputPGM p)
{
  scalar s = p.s;
  if (p.width == 0.) p.width = L0;

  char line[81];
  if (!fgets (line, 81, p.fp)) {
    fprintf (stderr, "input_pgm: could not read magic number\n");
    exit (1);
  }
  if (strcmp (line, "P5\n")) {
    fprintf (stderr, "input_pgm: magic number '%s' does not match PGM\n", 
	     line);
    exit (1);
  }
  if (!fgets (line, 81, p.fp)) {
    fprintf (stderr, "input_pgm: could not read width and height\n");
    exit (1);
  }
  int width, height;
  while (line[0] == '#' && fgets (line, 81, p.fp));
  if (line[0] == '#' || sscanf (line, "%d %d", &width, &height) != 2) {
    fprintf (stderr, "input_pgm: could not read width and height\n");
    exit (1);
  }
  if (!fgets (line, 81, p.fp)) {
    fprintf (stderr, "input_pgm: could not read maxval\n");
    exit (1);
  }
  int maxval;
  if (sscanf (line, "%d", &maxval) != 1) {
    fprintf (stderr, "input_pgm: could not read maxval\n");
    exit (1);
  }
  if (maxval < 256) {
    unsigned char * a = malloc (width*height);
    size_t n = fread (a, 1, width*height, p.fp);
    if (n != width*height) {
      fprintf (stderr, "input_pgm: read only %ld values\n", n);
      exit (1);
    }
    foreach() {
      int i = (x - p.ox)*width/p.width, j = (y - p.oy)*width/p.width;
      if (i >= 0 && i < width && j >= 0 && j < height)
	s[] = 1. - a[(height - 1 - j)*width + i]/(double)maxval;
      else
	s[] = 0.;
    }
    free (a);
  }
  else {
    unsigned short * a = malloc (width*height);
    size_t n = fread (a, 2, width*height, p.fp);
    if (n != width*height) {
      fprintf (stderr, "input_pgm: read only %ld values\n", n);
      exit (1);
    }
    foreach() {
      int i = (x - p.ox)*width/p.width, j = (y - p.oy)*width/p.width;
      if (i >= 0 && i < width && j >= 0 && j < height)
	s[] = 1. - a[(height - 1 - j)*width + i]/(double)maxval;
      else
	s[] = 0.;
    }
    free (a);
  }
}
