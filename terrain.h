#include <kdt/kdt.h>
#ifdef _OPENMP
# define NPROC omp_get_max_threads()
#else
# define NPROC 1
#endif

typedef struct {
  Kdt ** kdt;
  scalar n, dmin, dmax;
} Terrain;

Terrain * _terrain = NULL;
int _nterrain = 0;

static int includes (KdtRect rect, Point * p)
{
  Point point = *p;
  delta /= 2.;
  return (rect[0].l >= x - delta && rect[0].h <= x + delta &&
	  rect[1].l >= y - delta && rect[1].h <= y + delta);
}

static int intersects (KdtRect rect, Point * p)
{
  Point point = *p;
  delta /= 2.;
  return (rect[0].l <= x + delta && rect[0].h >= x - delta &&
	  rect[1].l <= y + delta && rect[1].h >= y - delta);
}

static void reconstruct_terrain (Point point, scalar zb)
{
  KdtSum s;
  kdt_sum_init (&s);
  delta /= 2.;
  KdtRect rect = {{x - delta, x + delta},
		  {y - delta, y + delta}};
  kdt_query_sum (_terrain[zb].kdt[pid()],
		 (KdtCheck) includes,
		 (KdtCheck) intersects, &point,
		 rect, &s);
  val(_terrain[zb].n,0,0) = s.n;
  if (s.w > 0.) {
    zb[] = s.H0/s.w;
    val(_terrain[zb].dmin,0,0) = s.Hmin;
    val(_terrain[zb].dmax,0,0) = s.Hmax;
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
    val(_terrain[zb].dmin,0,0) = nodata;
    val(_terrain[zb].dmax,0,0) = nodata;
  }
}

static void refine_terrain (Point point, scalar zb)
{
  foreach_child()
    reconstruct_terrain (point, zb);
}

void terrain (scalar zb, const char * name)
{
  if (zb >= _nterrain) {
    _terrain = realloc (_terrain, sizeof (Terrain)*(zb + 1));
    _nterrain = zb + 1;
  }
  _terrain[zb].kdt = malloc (NPROC*sizeof (Kdt *));
  for (int i = 0; i < NPROC; i++) {
    _terrain[zb].kdt[i] = kdt_new();
    if (kdt_open (_terrain[zb].kdt[i], name)) {
      fprintf (stderr, "terrain: could not open terrain database '%s'\n", name);
      exit (1);
    }
  }

  zb.refine = refine_terrain;
  scalar n = new scalar;
  n.refine = refine_none;
  n.coarsen = refine_none;
  _terrain[zb].n = n;
  scalar dmin = new scalar;
  dmin.refine = refine_none;
  dmin.coarsen = refine_none;
  _terrain[zb].dmin = dmin;
  scalar dmax = new scalar;
  dmax.refine = refine_none;
  dmax.coarsen = refine_none;
  _terrain[zb].dmax = dmax;

  trash ({zb, n, dmin, dmax});
  for (int l = 0; l <= depth(); l++) {
    foreach_level (l)
      reconstruct_terrain (point, zb);
    boundary_level ({zb}, l);
  }
}
