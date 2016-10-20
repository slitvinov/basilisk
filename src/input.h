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

#if MULTIGRID && !MULTIGRID_MPI

static void next_char (FILE * fp, int target)
{
  int c = fgetc(fp), para = 0;
  while (c != EOF && (c != target || para > 0)) {
    if (c == '{') para++;
    if (c == '}') para--;
    c = fgetc(fp);
  }
  if (c != target) {
    fprintf (stderr, "input_gfs(): error: expecting '%c'\n", target);
    exit (1);
  }
}

static int next_string (FILE * fp, const char * target)
{
  int slen = strlen (target), para = 0;
  char * s = malloc (slen + 1);
  s[slen] = '\0';
  int len = 0, c = fgetc (fp);
  while (c != EOF && len < slen) {
    if (c == '{') para++;
    if (c == '}') para--;
    s[len++] = c;
    c = fgetc (fp);
  }
  while (c != EOF && para >= 0) {
    if (!strcmp (s, target) && para == 0)
      break;
    if (c == '{') para++;
    if (c == '}') para--;
    for (int i = 0; i < slen - 1; i++)
      s[i] = s[i+1];
    s[slen - 1] = c;
    c = fgetc (fp);
  }
  if (strcmp (s, target))
    c = -1;
  free (s);
  return c;
}

/**
# *input_gfs()*: Gerris simulation format

The function reads simulation data in the format used in
[Gerris](http://gfs.sf.net) simulation files. This is the reciprocal
function of [*output_gfs()*](output.h#...).

The arguments and their default values are:

*fp*
: a file pointer. Default is *file* or stdin.

*list*
: a list of scalar fields to read. Default is *all*.

*file*
: the name of the file to read from (mutually exclusive with *fp*). */

void input_gfs (struct OutputGfs p)
{
  bool opened = false;
  if (p.fp == NULL) {
    if (p.file == NULL)
      p.fp = stdin;
    else if (!(p.fp = fopen (p.file, "r"))) {
      perror (p.file);
      exit (1);
    }
    else
      opened = true;
  }
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
    char * name = replace (s1, '_', '.', false);
    bool found = false;
    for (scalar s in p.list) {
      if (!found && !strcmp (s.name, name)) {
	input = list_append (input, s);
	found = true;
      }
    }
    if (!found)
      input = list_append (input, (scalar){INT_MAX});
    free (name);
    s1 = strtok (NULL, ", \t");
  }
  free (s);

  next_char (p.fp, '{');
  double t1 = 0.;
  if (next_string (p.fp, "Time") >= 0) {
    next_char (p.fp, '{');
    next_char (p.fp, 't');
    next_char (p.fp, '=');
    if (fscanf (p.fp, "%lf", &t1) != 1) {
      fprintf (stderr, "input_gfs(): error: expecting 't'\n");
      exit (1);
    }
    next_char (p.fp, '}');
    next_char (p.fp, '}');
  }

  if (next_string (p.fp, "Box") < 0) {
    fprintf (stderr, "input_gfs(): error: expecting 'GfsBox'\n");
    exit (1);
  }

  next_char (p.fp, '{');
  next_char (p.fp, '{');
  next_char (p.fp, '\n');

  scalar * listm = {cm,fm};
  scalar * listr = !is_constant(cm) ? listm : NULL;
  NOT_UNUSED (listr);

  #if TREE
  init_grid (1);
  #endif
  foreach_cell() {
    unsigned flags;
    if (fread (&flags, sizeof (unsigned), 1, p.fp) != 1) {
      fprintf (stderr, "input_gfs(): error: expecting 'flags'\n");
      exit (1);
    }
    if (!(flags & (1 << 4)) && is_leaf(cell))
      refine_cell (point, listr, 0, NULL);
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
      if (s.i != INT_MAX) {
	if (s.v.x.i >= 0) {
	  // this is a vector component, we need to rotate from
	  // Z-ordering (Gerris) to N-ordering (Basilisk)
#if dimension >= 2
	  if (s.v.x.i == s.i) {
	    s = s.v.y;
	    s[] = a;
	  }
	  else if (s.v.y.i == s.i) {
	    s = s.v.x;
	    s[] = - a;
	  }
#endif
#if dimension >= 3
	  else
	    s[] = a;
#endif
	}
	else	
	  s[] = a;
      }
    }
    if (is_leaf(cell))
      continue;
  }
  boundary (listm);
  boundary (input);

  free (input);
  if (opened)
    fclose (p.fp);

  // the events are advanced to catch up with the time
  while (t < t1 && events (false))
    t = tnext;
  events (false);
}

#endif // MULTIGRID && !MULTIGRID_MPI

/**
# *input_grd()*: Raster format (Esri grid)

This function reads a scalar field from a [Raster
file](https://en.wikipedia.org/wiki/Esri_grid). This is the reciprocal
function of [*output_grd()*](output.h#...).

The arguments and their default values are:

*s* 
: the scalar where the data will be stored. No default value. You
must specify this parameter

*fp*
: a file pointer. Default is *file* or stdin.

*file*
: the name of the file to read from (mutually exclusive with *fp*).

*nodatavalue*
: the value of the NoDataValue. Default is the same as that defined in
the raster file. 

*linear*
: if true, the raster data is bilinearly interpolated. Default is false.
*/

struct InputGRD {
  scalar s;
  FILE * fp;
  char * file;
  double nodatavalue;
  bool linear;
};

void input_grd (struct InputGRD p)
{
  scalar input = p.s;

  bool opened = false;
  if (p.fp == NULL) {
    if (p.file == NULL)
      p.fp = stdin;
    else if (!(p.fp = fopen (p.file, "r"))) {
      perror (p.file);
      exit (1);
    }
    else
      opened = true;
  }
  
  // Variables for the Raster data
  double DeltaGRD;
  int nx, ny;
  double XG0, YG0, ndv;

  // header
  char waste[100];
  fscanf (p.fp, "%s %d", waste, &nx);
  fscanf (p.fp, "%s %d", waste, &ny);
  fscanf (p.fp, "%s %lf", waste, &XG0);
  fscanf (p.fp, "%s %lf", waste, &YG0);
  fscanf (p.fp, "%s %lf", waste, &DeltaGRD);
  fscanf (p.fp, "%s %lf", waste, &ndv);

  //default value of NoData value
  if (!p.nodatavalue)
    p.nodatavalue = ndv;

  // read the data
  double * value = malloc (sizeof(double)*nx*ny);
  for (int i = ny - 1; i >= 0; i--)
    for (int j = 0 ; j < nx; j++)
      fscanf (p.fp, "%lf ", &value[j + i*nx]);

  double LGx0 = nx*DeltaGRD;
  double LGy0 = ny*DeltaGRD;
  bool warning = false;
  
  foreach() {
    // Test if the point in the Basilisk area is included in the raster area
    bool incl = (x >= XG0 - DeltaGRD/2. &&
		 x <= XG0 + LGx0 + DeltaGRD/2. &&
		 y >= YG0 - DeltaGRD/2. &&
		 y <= YG0 + LGy0 + DeltaGRD/2.);
    if (incl) {
      double val;
      // Test if we are on the ring of data around the raster grid
      bool ring = (x >= XG0 &&
		   x <= XG0 + LGx0 &&
		   y >= YG0 &&
		   y <= YG0 + LGy0);      
      if (p.linear && ring ) { // bi-linear interpolation
	int j = (x - XG0)/DeltaGRD; int i = (y - YG0)/DeltaGRD;
	double dx = x -(j*DeltaGRD + XG0); double dy = y - (i*DeltaGRD + YG0);
	val = value[j + i*nx]
	  + dx*(value[j + 1 + i*nx] - value[j + i*nx])/DeltaGRD
	  + dy*(value[j + (i + 1)*nx] - value[j + i*nx])/DeltaGRD
	  + dx*dy*(value[j + i*nx] + value[j +1 + (i+1)*nx]
		   - value[j + (i + 1)*nx]-value[j + 1 + i*nx])/sq(DeltaGRD);
      }
      else { 
	int j = (x - XG0 + DeltaGRD/2.)/DeltaGRD;
	int i = (y - YG0 + DeltaGRD/2.)/DeltaGRD;
	val = value[j + i*nx];
      }
      input[] = val;
      if (val == ndv)
	input[] = p.nodatavalue;
    }
    else {
      input[] = p.nodatavalue;
      warning = true;
    }
  }
  
  if (warning)
    fprintf (stderr,
	     "input_grd(): Warning: Raster data is not covering all"
	     " the simulation area\n");

  free (value);
  if (opened)
    fclose (p.fp);
}
