struct OutputField {
  scalar * list;
  FILE * fp;
  int n;
  bool linear;
};

void output_field (struct OutputField p)
{
  fprintf (p.fp, "# 1:x 2:y");
  int i = 3;
  for (scalar s in p.list)
    fprintf (p.fp, " %d:%d", i++, s);
  fputc('\n', p.fp);
  double delta = L0/p.n;
  for (int i = 0; i < p.n; i++) {
    double xp = delta*i + X0 + delta/2.;
    for (int j = 0; j < p.n; j++) {
      double yp = delta*j + Y0 + delta/2.;
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

struct OutputMatrix {
  scalar f;
  FILE * fp;
  int n;
  bool linear;
};

void output_matrix (struct OutputMatrix p)
{
  float fn = p.n;
  float delta = L0/fn;
  fwrite (&fn, sizeof(float), 1, p.fp);
  for (int j = 0; j < p.n; j++) {
    float yp = delta*j + X0 + delta/2.;
    fwrite (&yp, sizeof(float), 1, p.fp);
  }
  for (int i = 0; i < p.n; i++) {
    float xp = delta*i + X0 + delta/2.;
    fwrite (&xp, sizeof(float), 1, p.fp);
    for (int j = 0; j < p.n; j++) {
      float yp = delta*j + Y0 + delta/2., v;
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

struct OutputPPM {
  scalar f;
  FILE * fp;
  int n;
  double min, max;
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
    double avg = s.sum/s.area;
    p.min = avg - 10.*s.stddev; p.max = avg + 10.*s.stddev;
  }
  if (p.box[0][0] == 0. && p.box[0][1] == 0. && 
      p.box[1][0] == 0. && p.box[1][1] == 0.) {
    p.box[0][0] = X0;      p.box[0][1] = Y0;
    p.box[1][0] = X0 + L0; p.box[1][1] = Y0 + L0;
  }

  double fn = p.n;
  double delta = (p.box[1][0] - p.box[0][0])/fn;
  int ny = (p.box[1][1] - p.box[0][1])/delta;
  
  color ** ppm = matrix_new (ny, p.n, sizeof(color));
  double cmap[NCMAP][3];
  colormap_jet (cmap);
  OMP_PARALLEL()
  OMP(omp for schedule(static))
  for (int j = 0; j < ny; j++) {
    double yp = delta*j + p.box[0][1] + delta/2.;
    for (int i = 0; i < p.n; i++) {
      double xp = delta*i + p.box[0][0] + delta/2., v;
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

  fprintf (p.fp, "P6\n%u %u 255\n", p.n, ny);
  fwrite (((void **) ppm)[0], sizeof(color), ny*p.n, p.fp);
  fflush (p.fp);
  
  matrix_free (ppm);
}