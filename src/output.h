/**
# Output functions

## Multiple fields interpolated on a regular grid (text format)

This function interpolates a *list* of fields on a *n x n* regular
grid. The resulting data are written in text format in the file
pointed to by *fp*. The correspondance between column numbers and
variables is summarised in the first line of the file. The data are
written row-by-row and each row is separated from the next by a
blank line. This format is compatible with the *splot* command of
*gnuplot* i.e. one could use something like

~~~bash
gnuplot> set pm3d map
gnuplot> splot 'fields' u 1:2:4
~~~

The arguments and their default values are:

*list*
: list of fields to output. Default is *all*.

*fp*
: file pointer. Default is *stdout*.

*n*
: number of points along each dimension. Default is *N*.

*linear*
: use first-order (default) or bilinear interpolation. */

struct OutputField {
  scalar * list;
  FILE * fp;
  int n;
  bool linear;
};

void output_field (struct OutputField p)
{
  if (!p.list) p.list = all;
  if (p.n == 0) p.n = N;
  if (!p.fp) p.fp = stdout;
  fprintf (p.fp, "# 1:x 2:y");
  int i = 3;
  for (scalar s in p.list)
    fprintf (p.fp, " %d:%s", i++, s.name);
  fputc('\n', p.fp);
  double Delta = L0/p.n;
  for (int i = 0; i < p.n; i++) {
    double xp = Delta*i + X0 + Delta/2.;
    for (int j = 0; j < p.n; j++) {
      double yp = Delta*j + Y0 + Delta/2.;
      fprintf (p.fp, "%g %g", xp, yp);
      if (p.linear) {
	for (scalar s in p.list)
	  fprintf (p.fp, " %g", interpolate (s, xp, yp));
      }
      else {
	Point point = locate (xp, yp);
	for (scalar s in p.list)
	  fprintf (p.fp, " %g", point.level >= 0 ? s[] : nodata);
      }
      fputc ('\n', p.fp);
    }
    fputc ('\n', p.fp);
  }
  fflush (p.fp);
}

/**
## Single field interpolated on a regular grid (binary format)

This function writes a binary representation of a single field
interpolated on a regular *n x n* grid. The format is compatible with
the binary matrix format of gnuplot i.e. one could use

~~~bash
gnuplot> set pm3d map
gnuplot> splot 'matrix' binary u 2:1:3
~~~

The arguments and their default values are:

*f*
: a scalar field (compulsory).

*fp*
: file pointer. Default is *stdout*.

*n*
: number of points along each dimension. Default is *N*.

*linear*
: use first-order (default) or bilinear interpolation. */

struct OutputMatrix {
  scalar f;
  FILE * fp;
  int n;
  bool linear;
};

void output_matrix (struct OutputMatrix p)
{
  if (p.n == 0) p.n = N;
  if (!p.fp) p.fp = stdout;
  float fn = p.n;
  float Delta = L0/fn;
  fwrite (&fn, sizeof(float), 1, p.fp);
  for (int j = 0; j < p.n; j++) {
    float yp = Delta*j + X0 + Delta/2.;
    fwrite (&yp, sizeof(float), 1, p.fp);
  }
  for (int i = 0; i < p.n; i++) {
    float xp = Delta*i + X0 + Delta/2.;
    fwrite (&xp, sizeof(float), 1, p.fp);
    for (int j = 0; j < p.n; j++) {
      float yp = Delta*j + Y0 + Delta/2., v;
      if (p.linear)
	v = interpolate (p.f, xp, yp);
      else {
	Point point = locate (xp, yp);
	assert (point.level >= 0);
	v = val(p.f,0,0);
      }
      fwrite (&v, sizeof(float), 1, p.fp);
    }
  }
  fflush (p.fp);
}

/**
## Colormaps

Colormaps are arrays of (127) red, green, blue triplets. For the
moment we only use the "jet" colormap. */

#define NCMAP 127

void colormap_jet (double cmap[NCMAP][3])
{
  for (int i = 0; i < NCMAP; i++) {
    cmap[i][0] = 
      i <= 46 ? 0. : 
      i >= 111 ? -0.03125*(i - 111) + 1. :
      i >= 78 ? 1. : 
      0.03125*(i - 46);
    cmap[i][1] = 
      i <= 14 || i >= 111 ? 0. : 
      i >= 79 ? -0.03125*(i - 111) : 
      i <= 46 ? 0.03125*(i - 14) : 
      1.;
    cmap[i][2] =
      i >= 79 ? 0. :
      i >= 47 ? -0.03125*(i - 79) :
      i <= 14 ? 0.03125*(i - 14) + 1.:
      1.;
  }
}

/**
Given a colormap, a minimum and maximum value, the this function
returns the red/green/blue triplet corresponding to *val*. */

typedef struct {
  unsigned char r, g, b;
} color;

color colormap_color (double cmap[NCMAP][3], 
		      double val, double min, double max)
{
  color c;
  if (val == nodata) {
    c.r = c.g = c.b = 0; // nodata is black
    return c;
  }    
  val = val <= min ? 0. : val >= max ? 0.9999 : (val - min)/(max - min);
  int i = val*(NCMAP - 1);
  double coef = val*(NCMAP - 1) - i;
  assert (i < NCMAP - 1);
  unsigned char * c1 = (unsigned char *) &c;
  for (int j = 0; j < 3; j++)
    c1[j] = 255*(cmap[i][j]*(1. - coef) + cmap[i + 1][j]*coef);
  return c;
}

/**
## Portable PixMap (PPM) image output

Given a field, this function outputs a colormaped representation as a
[Portable PixMap](http://en.wikipedia.org/wiki/Netpbm_format) image.

If [ImageMagick](http://www.imagemagick.org/) is installed on the
system, this image can optionally be converted to any image format
supported by ImageMagick.

The arguments and their default values are:

*f*
: a scalar field (compulsory).

*fp*
: a file pointer. Default is stdout.

*n*
: number of pixels. Default is *N*.

*file*
: sets the name of the file used as output for
ImageMagick. This allows outputs in all formats supported by
ImageMagick. For example, one could use

~~~c
output_ppm (f, file = "f.png");
~~~

to get a [PNG](http://en.wikipedia.org/wiki/Portable_Network_Graphics)
image.

*min, max*
: minimum and maximum values used to define the
colorscale. By default these are set automatically using the *spread*
parameter. 

*spread*
: if not specified explicitly, *min* and *max* are set to the average 
of the field minus (resp. plus) *spread* times the standard deviation. 
By default *spread* is five. 

*linear*
: whether to use bilinear or first-order interpolation. Default is 
first-order.

*box*
: the lower-left and upper-right coordinates of the domain to consider.
 Default is the entire domain.

*mask*
: if set, this field will be used to mask out (in black), the regions 
of the domain for which *mask* is negative. */

struct OutputPPM {
  scalar f;
  FILE * fp;
  int n;
  char * file;
  double min, max, spread;
  bool linear;
  double box[2][2];
  scalar mask;
};

void output_ppm (struct OutputPPM p)
{
  // default values
  if (p.n == 0) p.n = N;
  if (p.min == 0 && p.max == 0) {
    stats s = statsf (p.f);
    double avg = s.sum/s.area, spread = (p.spread ? p.spread : 5.)*s.stddev;
    p.min = avg - spread; p.max = avg + spread;
  }
  if (p.box[0][0] == 0. && p.box[0][1] == 0. && 
      p.box[1][0] == 0. && p.box[1][1] == 0.) {
    p.box[0][0] = X0;      p.box[0][1] = Y0;
    p.box[1][0] = X0 + L0; p.box[1][1] = Y0 + L0;
  }

  double fn = p.n;
  double Delta = (p.box[1][0] - p.box[0][0])/fn;
  int ny = (p.box[1][1] - p.box[0][1])/Delta;
  
  color ** ppm = matrix_new (ny, p.n, sizeof(color));
  double cmap[NCMAP][3];
  colormap_jet (cmap);
  OMP_PARALLEL()
  OMP(omp for schedule(static))
  for (int j = 0; j < ny; j++) {
    double yp = Delta*j + p.box[0][1] + Delta/2.;
    for (int i = 0; i < p.n; i++) {
      double xp = Delta*i + p.box[0][0] + Delta/2., v;
      if (p.mask) { // masking
	if (p.linear) {
	  double m = interpolate (p.mask, xp, yp);
	  if (m < 0.)
	    v = nodata;
	  else
	    v = interpolate (p.f, xp, yp);
	}
	else {
	  Point point = locate (xp, yp);
	  if (point.level < 0 || val(p.mask,0,0) < 0.)
	    v = nodata;
	  else
	    v = val(p.f,0,0);
	}
      }
      else if (p.linear)
	v = interpolate (p.f, xp, yp);
      else {
	Point point = locate (xp, yp);
	v = point.level >= 0 ? val(p.f,0,0) : nodata;
      }
      ppm[ny - 1 - j][i] = colormap_color (cmap, v, p.min, p.max);
    }
  }
  OMP_END_PARALLEL()

  if (!p.fp) p.fp = stdout;
  if (p.n == 0) p.n = N;
  if (p.file) {
    char * command = malloc (strlen ("convert ppm:- ") + strlen (p.file) + 1);
    strcpy (command, "convert ppm:- ");
    strcat (command, p.file);
    p.fp = popen (command, "w");
    free (command);
  }

  fprintf (p.fp, "P6\n%u %u 255\n", p.n, ny);
  fwrite (((void **) ppm)[0], sizeof(color), ny*p.n, p.fp);

  if (p.file)
    pclose (p.fp);
  else
    fflush (p.fp);
  
  matrix_free (ppm);
}

/**
## ESRI ASCII Grid format

The [ESRI GRD format](http://en.wikipedia.org/wiki/Esri_grid) is a
standard format for importing raster data into [GIS
systems](http://en.wikipedia.org/wiki/Geographic_information_system).

The arguments and their default values are:

*f*
: a scalar field (compulsory).

*fp*
: a file pointer. Default is stdout.

$\Delta$
: size of a grid element. Default is 1/*N*.

*linear*
: whether to use bilinear or first-order interpolation. Default is 
first-order.

*box*
: the lower-left and upper-right coordinates of the domain to consider.
 Default is the entire domain.

*mask*
: if set, this field will be used to mask out, the regions 
of the domain for which *mask* is negative. */

struct OutputGRD {
  scalar f;
  FILE * fp;
  double Delta;
  bool linear;
  double box[2][2];
  scalar mask;
};

void output_grd (struct OutputGRD p)
{
  // default values
  if (!p.fp) p.fp = stdout;
  if (p.box[0][0] == 0. && p.box[0][1] == 0. && 
      p.box[1][0] == 0. && p.box[1][1] == 0.) {
    p.box[0][0] = X0;      p.box[0][1] = Y0;
    p.box[1][0] = X0 + L0; p.box[1][1] = Y0 + L0;
    if (p.Delta == 0) p.Delta = L0/N;
  }

  double Delta = p.Delta;
  int nx = (p.box[1][0] - p.box[0][0])/Delta;
  int ny = (p.box[1][1] - p.box[0][1])/Delta;

  // header
  fprintf (p.fp, "ncols          %d\n", nx);
  fprintf (p.fp, "nrows          %d\n", ny);
  fprintf (p.fp, "xllcorner      %g\n", p.box[0][0]);
  fprintf (p.fp, "yllcorner      %g\n", p.box[0][1]);
  fprintf (p.fp, "cellsize       %g\n", Delta);
  fprintf (p.fp, "nodata_value   -9999\n");
  
  // data
  for (int j = ny-1; j >= 0; j--) {
    double yp = Delta*j + p.box[0][1] + Delta/2.;
    for (int i = 0; i < nx; i++) {
      double xp = Delta*i + p.box[0][0] + Delta/2., v;
      if (p.mask) { // masking
	if (p.linear) {
	  double m = interpolate (p.mask, xp, yp);
	  if (m < 0.)
	    v = nodata;
	  else
	    v = interpolate (p.f, xp, yp);
	}
	else {
	  Point point = locate (xp, yp);
	  if (point.level < 0 || val(p.mask,0,0) < 0.)
	    v = nodata;
	  else
	    v = val(p.f,0,0);
	}
      }
      else if (p.linear)
	v = interpolate (p.f, xp, yp);
      else {
	Point point = locate (xp, yp);
	v = point.level >= 0 ? val(p.f,0,0) : nodata;
      }
      if (v == nodata)
	fprintf (p.fp, "-9999 ");
      else
	fprintf (p.fp, "%f ", v);
    }
    fprintf (p.fp, "\n");
  }

  fflush (p.fp);
}
