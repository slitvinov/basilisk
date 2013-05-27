void output_field (scalar * list, int n, FILE * fp, bool linear)
{
  fprintf (fp, "# 1:x 2:y");
  int i = 3;
  for (scalar s in list)
    fprintf (fp, " %d:%d", i++, s);
  fputc('\n', fp);
  double delta = L0/n;
  for (int i = 0; i < n; i++) {
    double x = delta*i + X0 + delta/2.;
    for (int j = 0; j < n; j++) {
      double y = delta*j + Y0 + delta/2.;
      fprintf (fp, "%g %g", x, y);
      if (linear) {
	for (scalar s in list)
	  fprintf (fp, " %g", interpolate (s, x, y));
      }
      else {
	Point point = locate (x, y);
	for (scalar s in list)
	  fprintf (fp, " %g", point.level >= 0 ? val(s,0,0) : nodata);
      }
      fputc ('\n', fp);
    }
    fputc ('\n', fp);
  }
  fflush (fp);
}

void output_matrix (scalar f, int n, FILE * fp, bool linear)
{
  float fn = n;
  float delta = L0/fn;
  fwrite (&fn, sizeof(float), 1, fp);
  for (int j = 0; j < n; j++) {
    float y = delta*j + X0 + delta/2.;
    fwrite (&y, sizeof(float), 1, fp);
  }
  for (int i = 0; i < n; i++) {
    float x = delta*i + X0 + delta/2.;
    fwrite (&x, sizeof(float), 1, fp);
    for (int j = 0; j < n; j++) {
      float y = delta*j + Y0 + delta/2., v;
      if (linear)
	v = interpolate (f, x, y);
      else {
	Point point = locate (x, y);
	assert (point.level >= 0);
	v = val(f,0,0);
      }
      fwrite (&v, sizeof(float), 1, fp);
    }
  }
  fflush (fp);
}

#define NCMAP 127

void colormap_jet (float cmap[NCMAP][3])
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

void colormap_color (float cmap[NCMAP][3], 
		     float val, float min, float max,
		     unsigned char c[3])
{
  val = val <= min ? 0. : val >= max ? 0.9999 : (val - min)/(max - min);
  int i = val*(NCMAP - 1);
  float coef = val*(NCMAP - 1) - i;
  assert (i < NCMAP - 1);
  for (int j = 0; j < 3; j++)
    c[j] = 255*(cmap[i][j]*(1. - coef) + cmap[i + 1][j]*coef);
}

struct OutputPPM {
  scalar f;
  FILE * fp;
  int n;
  double min, max;
  bool linear;
  double box[2][2];
};

void output_ppm (struct OutputPPM p)
{
  // default values
  if (p.n == 0) p.n = N;
  if (p.min == 0 && p.max == 0) {
    stats s = statsf (p.f);
    p.min = s.min; p.max = s.max;
  }
  if (p.box[0][0] == 0. && p.box[0][1] == 0. && 
      p.box[1][0] == 0. && p.box[1][1] == 0.) {
    p.box[0][0] = X0;      p.box[0][1] = Y0;
    p.box[1][0] = X0 + L0; p.box[1][1] = Y0 + L0;
  }

  double fn = p.n;
  double delta = (p.box[1][0] - p.box[0][0])/fn;
  int ny = (p.box[1][1] - p.box[0][1])/delta;
  fprintf (p.fp, "P6\n%u %u 255\n", p.n, ny);
  float cmap[NCMAP][3];
  colormap_jet (cmap);
  for (int j = ny - 1; j >= 0; j--) {
    double y = delta*j + p.box[0][1] + delta/2.;
    for (int i = 0; i < p.n; i++) {
      double x = delta*i + p.box[0][0] + delta/2., v;
      if (p.linear)
	v = interpolate (p.f, x, y);
      else {
	Point point = locate (x, y);
	v = point.level >= 0 ? val(p.f,0,0) : nodata;
      }
      unsigned char c[3];
      if (v == nodata)
	c[0] = c[1] = c[2] = 0;
      else
	colormap_color (cmap, v, p.min, p.max, c);
      fwrite (c, sizeof(unsigned char), 3, p.fp);
    }
  }
  fflush (p.fp);
}
