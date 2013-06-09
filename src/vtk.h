void output_vtk (scalar * list, int n, FILE * fp, bool linear)
{
  fputs ("# vtk DataFile Version 2.0\n"
	 "Basilisk\n"
	 "ASCII\n"
	 "DATASET STRUCTURED_GRID\n", fp);
  fprintf (fp, "DIMENSIONS %d %d 1\n", n, n);
  fprintf (fp, "POINTS %d double\n", n*n);

  double fn = n;
  double delta = L0/fn;
  for (int i = 0; i < n; i++) {
    double x = delta*i + X0 + delta/2.;
    for (int j = 0; j < n; j++) {
      double y = delta*j + Y0 + delta/2.;
      fprintf (fp, "%g %g 0\n", x, y);
    }
  }
  fprintf (fp, "POINT_DATA %d\n", n*n);
  for (scalar s in list) {
    fprintf (fp, "SCALARS scalar%d double\n", s);
    fputs ("LOOKUP_TABLE default\n", fp);
    double fn = n;
    double delta = L0/fn;
    for (int i = 0; i < n; i++) {
      double x = delta*i + X0 + delta/2.;
      for (int j = 0; j < n; j++) {
	double y = delta*j + Y0 + delta/2., v;
	if (linear)
	  v = interpolate (s, x, y);
	else {
	  Point point = locate (x, y);
	  v = point.level >= 0 ? val(s,0,0) : nodata;
	}
	fprintf (fp, "%g\n", v);
      }
    }
  }
  fflush (fp);
}
