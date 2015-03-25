/* tangential interpolation on face vector fields (in parallel) 
 see also: interpu.c */

scalar h[];
face vector u[];

double R0 = 0.1;
u.n[right] = dirichlet(exp(-(x*x + y*y)/(R0*R0)));
u.t[top] = dirichlet(exp(-(x*x + y*y)/(R0*R0)));

int main (int argc, char ** argv)
{
  int maxlevel = argc > 1 ? atoi(argv[1]) : 6;
  int minlevel = argc > 2 ? atoi(argv[2]) : 5;
  origin (-1, -1);
  init_grid (1);
  refine (level <= minlevel*(1. - sqrt(sq((x + 0.5) - 0.1) +
				       sq((y + 0.5)- 0.1))), NULL);
  mpi_partitioning();
  
  refine (level <= maxlevel*(1. - sqrt(sq((x + 0.5) - 0.1) +
				       sq((y + 0.5) - 0.1))), NULL);
  
  foreach()
    h[] = exp(-(x*x + y*y)/(R0*R0));
  boundary ({h});
  
  foreach_face(x) u.x[] = exp(-(x*x + y*y)/(R0*R0));
  foreach_face(y) u.y[] = exp(-(x*x + y*y)/(R0*R0));
  boundary ((scalar *){u});
  //  output_cells (stdout);

  double max = 0;
  foreach_face (x)
    for (int i = -1; i <= 1; i++) {
      double xu = x, yu = y + i*Delta;
      double e = exp(-(xu*xu+yu*yu)/(R0*R0)) - u.x[0,i];
      if (fabs(e) > max)
	max = fabs(e);
    }

  double maxv = 0;
  foreach_face (y)
    for (int i = -1; i <= 1; i++) {
      double xv = x + i*Delta, yv = y;
      double e = exp(-(xv*xv+yv*yv)/(R0*R0)) - u.y[i,0];
      if (fabs(e) > maxv)
	maxv = fabs(e);
      printf ("%g %g %d %d %g %g\n", 
	      xv, yv, level, cell.neighbors, u.y[i,0], e);
    }

  fprintf (stderr, "maximum error on halos: %g %g\n", max, maxv);
}
