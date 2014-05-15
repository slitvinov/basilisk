#include <stdarg.h>
#include <kdt/kdt.h>
#ifdef _OPENMP
# define NPROC omp_get_max_threads()
#else
# define NPROC 1
#endif

attribute {
  void ** kdt;
  scalar n, dmin, dmax;
};

static int includes (KdtRect rect, Point * p)
{
  Point point = *p;
  Delta /= 2.;
  return (rect[0].l >= x - Delta && rect[0].h <= x + Delta &&
	  rect[1].l >= y - Delta && rect[1].h <= y + Delta);
}

static int intersects (KdtRect rect, Point * p)
{
  Point point = *p;
  Delta /= 2.;
  return (rect[0].l <= x + Delta && rect[0].h >= x - Delta &&
	  rect[1].l <= y + Delta && rect[1].h >= y - Delta);
}

static void reconstruct_terrain (Point point, scalar zb)
{
  KdtSum s;
  kdt_sum_init (&s);
  Delta /= 2.;
  KdtRect rect = {{x - Delta, x + Delta},
		  {y - Delta, y + Delta}};
  for (Kdt ** kdt = zb.kdt[pid()]; *kdt; kdt++)
    kdt_query_sum (*kdt,
		   (KdtCheck) includes,
		   (KdtCheck) intersects, &point,
		   rect, &s);
  val(zb.n,0,0) = s.n;
  if (s.w > 0.) {
    zb[] = s.H0/s.w;
    val(zb.dmin,0,0) = s.Hmin;
    val(zb.dmax,0,0) = s.Hmax;
  }
  else {
    /* not enough points in database, use bilinear interpolation
       from coarser level instead */
    if (level > 0)
      zb[] = (9.*coarse(zb,0,0) + 
	      3.*(coarse(zb,child.x,0) + coarse(zb,0,child.y)) + 
	      coarse(zb,child.x,child.y))/16.;
    else
      zb[] = 0.; // no points at level 0!
    val(zb.dmin,0,0) = nodata;
    val(zb.dmax,0,0) = nodata;
  }
}

static void refine_terrain (Point point, scalar zb)
{
  foreach_child()
    reconstruct_terrain (point, zb);
}

void terrain (scalar zb, ...)
{
  zb.kdt = calloc (NPROC, sizeof (Kdt **));

  int nt = 0;
  va_list ap;
  va_start (ap, zb);
  char * name;
  while ((name = va_arg (ap, char *))) {
    for (int i = 0; i < NPROC; i++) {
      Kdt ** kdt = zb.kdt[i];
      zb.kdt[i] = kdt = realloc (kdt, sizeof(Kdt *)*(nt + 2));
      kdt[nt] = kdt_new();
      kdt[nt + 1] = NULL;
      if (kdt_open (kdt[nt], name)) {
	fprintf (stderr, "terrain: could not open terrain database '%s'\n", 
		 name);
	exit (1);
      }
    }
    nt++;
  }
  va_end (ap);

  zb.refine = refine_terrain;
  scalar n = new scalar;
  n.refine = none;
  n.coarsen = none;
  zb.n = n;
  scalar dmin = new scalar;
  dmin.refine = none;
  dmin.coarsen = none;
  zb.dmin = dmin;
  scalar dmax = new scalar;
  dmax.refine = none;
  dmax.coarsen = none;
  zb.dmax = dmax;

  trash ({zb});
  for (int l = 0; l <= depth(); l++) {
    foreach_level (l)
      reconstruct_terrain (point, zb);
    boundary_level ({zb}, l);
  }
}
