/**
# Output functions

## *output_field()*: Multiple fields interpolated on a regular grid (text format)

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
## *output_matrix()*: Single field interpolated on a regular grid (binary format)

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

Colormaps are arrays of (127) red, green, blue triplets. */

#define NCMAP 127

typedef void (* colormap) (double cmap[NCMAP][3]);

void jet (double cmap[NCMAP][3])
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

void cool_warm (double cmap[NCMAP][3])
{
  /* diverging cool-warm from:
   *  http://www.sandia.gov/~kmorel/documents/ColorMaps/CoolWarmFloat33.csv
   * see also:
   *  Diverging Color Maps for Scientific Visualization (Expanded)
   *  Kenneth Moreland
   */
  static double basemap[33][3] = {
    {0.2298057,   0.298717966, 0.753683153},
    {0.26623388,  0.353094838, 0.801466763},
    {0.30386891,  0.406535296, 0.84495867},
    {0.342804478, 0.458757618, 0.883725899},
    {0.38301334,  0.50941904,  0.917387822},
    {0.424369608, 0.558148092, 0.945619588},
    {0.46666708,  0.604562568, 0.968154911},
    {0.509635204, 0.648280772, 0.98478814},
    {0.552953156, 0.688929332, 0.995375608},
    {0.596262162, 0.726149107, 0.999836203},
    {0.639176211, 0.759599947, 0.998151185},
    {0.681291281, 0.788964712, 0.990363227},
    {0.722193294, 0.813952739, 0.976574709},
    {0.761464949, 0.834302879, 0.956945269},
    {0.798691636, 0.849786142, 0.931688648},
    {0.833466556, 0.860207984, 0.901068838},
    {0.865395197, 0.86541021,  0.865395561},
    {0.897787179, 0.848937047, 0.820880546},
    {0.924127593, 0.827384882, 0.774508472},
    {0.944468518, 0.800927443, 0.726736146},
    {0.958852946, 0.769767752, 0.678007945},
    {0.96732803,  0.734132809, 0.628751763},
    {0.969954137, 0.694266682, 0.579375448},
    {0.966811177, 0.650421156, 0.530263762},
    {0.958003065, 0.602842431, 0.481775914},
    {0.943660866, 0.551750968, 0.434243684},
    {0.923944917, 0.49730856,  0.387970225},
    {0.89904617,  0.439559467, 0.343229596},
    {0.869186849, 0.378313092, 0.300267182},
    {0.834620542, 0.312874446, 0.259301199},
    {0.795631745, 0.24128379,  0.220525627},
    {0.752534934, 0.157246067, 0.184115123},
    {0.705673158, 0.01555616,  0.150232812}	
  };
  
  for (int i = 0; i < NCMAP; i++) {
    double x = i*(32 - 1e-10)/(NCMAP - 1);
    int j = x; x -= j;
    for (int k = 0; k < 3; k++)
      cmap[i][k] = (1. - x)*basemap[j][k] + x*basemap[j+1][k];
  }
}

void gray (double cmap[NCMAP][3])
{
  for (int i = 0; i < NCMAP; i++)
    for (int k = 0; k < 3; k++)
      cmap[i][k] = i/(NCMAP - 1.);
}

/**
Given a colormap and a minimum and maximum value, this function
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
## *output_ppm()*: Portable PixMap (PPM) image output

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
of the domain for which *mask* is negative. 

*map*
: the colormap: *jet*, *cool_warm* or *gray*. Default is *jet*.
*/

struct OutputPPM {
  scalar f;
  FILE * fp;
  int n;
  char * file;
  double min, max, spread;
  bool linear;
  double box[2][2];
  scalar mask;
  colormap map;
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
  if (!p.map)
    p.map = jet;

  double fn = p.n;
  double Delta = (p.box[1][0] - p.box[0][0])/fn;
  int ny = (p.box[1][1] - p.box[0][1])/Delta;
  
  color ** ppm = matrix_new (ny, p.n, sizeof(color));
  double cmap[NCMAP][3];
  p.map (cmap);
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
## *output_grd()*: ESRI ASCII Grid format

The [ESRI GRD format](http://en.wikipedia.org/wiki/Esri_grid) is a
standard format for importing raster data into [GIS
systems](http://en.wikipedia.org/wiki/Geographic_information_system).

The arguments and their default values are:

*f*
: a scalar field (compulsory).

*fp*
: a file pointer. Default is stdout.

$\Delta$
: size of a grid element. Default is 1/N.

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
