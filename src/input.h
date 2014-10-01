/**
# *input_pgm()*: Importing Portable Gray Map (PGM) images

This function reads in the
[PGM](http://en.wikipedia.org/wiki/Netpbm_format) file in file *fp*
and imports the corresponding values into field *s*. The grayscale is
normalised and inverted so that the maximum value in field *s* is one
(black) and the minimum value is zero (white). 

By default the origin of the image (lower-left corner) is assumed to
be at (0,0) and the width of the image is set to `L0`. This can be
changed using the optional *ox*, *oy* and *width* parameters. */

#include "utils.h"

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
  if (strcmp (line, "P2\n") && strcmp (line, "P5\n")) {
    fprintf (stderr, "input_pgm: magic number '%s' does not match PGM\n", 
	     line);
    exit (1);
  }
  int binary = !strcmp (line, "P5\n");
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
    size_t n = 0;
    if (binary)
      n = fread (a, 1, width*height, p.fp);
    else {
      int v;
      while (n < width*height && fscanf (p.fp, "%d ", &v) == 1)
	a[n++] = v;
    }
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
    unsigned short * a = malloc (2*width*height);
    size_t n = 0;
    if (binary)
      n = fread (a, 2, width*height, p.fp);
    else {
      int v;
      while (n < width*height && fscanf (p.fp, "%d ", &v) == 1)
	a[n++] = v;
    }
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

#if QUADTREE

static void next_char (FILE * fp, int target)
{
  int c = fgetc(fp);
  while (c != EOF && c != target)
    c = fgetc(fp);
  if (c != target) {
    fprintf (stderr, "input_gfs(): error: expecting '%c'\n", target);
    exit (1);
  }
}

/**
# *input_gfs()*: Gerris simulation format

The function reads simulation data in the format used in
[Gerris](http://gfs.sf.net) simulation files. This is the reciprocal
function of [*output_gfs()*](output.h#).

The arguments and their default values are:

*fp*
: a file pointer. Default is stdin.

*list*
: a list of scalar fields to read. Default is *all*. */

void input_gfs (struct OutputGfs p)
{
  if (p.fp == NULL) p.fp = stdin;
  if (p.list == NULL) p.list = all;

  next_char (p.fp, '{');
  
  char * s = malloc (1);
  int len = 0;
  int c = fgetc(p.fp);
  while (c != EOF && c != '}') {
    s[len++] = c;
    s = realloc (s, len + 1);
    s[len] = '\0';
    c = fgetc(p.fp);
  }
  if (c != '}') {
    fprintf (stderr, "input_gfs(): error: expecting '}'\n");
    exit (1);
  }
  
  char * s1 = strstr (s, "variables");
  if (!s1) {
    fprintf (stderr, "input_gfs(): error: expecting 'variables'\n");
    exit (1);
  }

  s1 = strstr (s1, "=");
  if (!s1) {
    fprintf (stderr, "input_gfs(): error: expecting '='\n");
    exit (1);
  }
  s1++;

  while (strchr (" \t", *s1))
    s1++;

  scalar * input = NULL;
  s1 = strtok (s1, ", \t");
  while (s1) {
    char * name = replace (s1, '_', '.');
    bool found = false;
    for (scalar s in p.list) {
      if (!found && !strcmp (s.name, name)) {
	input = list_append (input, s);
	found = true;
      }
    }
    if (!found)
      input = list_append (input, INT_MAX);
    free (name);
    s1 = strtok (NULL, ", \t");
  }

  s = realloc (s, 7);
  s[6] = '\0';
  len = 0;
  while (c != EOF && len < 6) {
    s[len++] = c;
    c = fgetc (p.fp);
  }
  while (c != EOF) {
    for (int i = 0; i < 5; i++)
      s[i] = s[i+1];
    s[5] = c;
    if (!strcmp (s, "GfsBox"))
      break;
    c = fgetc (p.fp);
  }
  if (strcmp (s, "GfsBox")) {
    fprintf (stderr, "input_gfs(): error: expecting 'GfsBox'\n");
    exit (1);
  }
  free (s);

  next_char (p.fp, '{');
  next_char (p.fp, '{');
  next_char (p.fp, '\n');

  scalar * listm = {cm,fm};
  scalar * listr = !is_constant(cm) ? listm : NULL;

  init_grid (1);
  foreach_cell() {
    unsigned flags;
    if (fread (&flags, sizeof (unsigned), 1, p.fp) != 1) {
      fprintf (stderr, "input_gfs(): error: expecting 'flags'\n");
      exit (1);
    }
    if (!(flags & (1 << 4)) && is_leaf(cell))
      point = refine_cell (point, listr, 0);
    double a;
    if (fread (&a, sizeof (double), 1, p.fp) != 1 || a != -1) {
      fprintf (stderr, "input_gfs(): error: expecting '-1'\n");
      exit (1);
    }
    for (scalar s in input) {
      if (fread (&a, sizeof (double), 1, p.fp) != 1) {
	fprintf (stderr, "input_gfs(): error: expecting a scalar\n");
	exit (1);
      }
      if (s != INT_MAX)
	s[] = a;
    }
    if (is_leaf(cell))
      continue;
  }
  boundary (listm);
  boundary (input);

  free (input);
}
#endif
