#include <kdt.h>

typedef struct {
  Kdt * kdt;
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

static void refine_terrain (Point point, scalar zb)
{
  foreach_child() {
    KdtSum s;
    kdt_sum_init (&s);
    delta /= 2.;
    KdtRect rect = {{x - delta, x + delta},
		    {y - delta, y + delta}};
    kdt_query_sum (_terrain[zb].kdt,
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
      zb[] = (9.*coarse(zb,0,0) + 
	      3.*(coarse(zb,child.x,0) + coarse(zb,0,child.y)) + 
	      coarse(zb,child.x,child.y))/16.;
      val(_terrain[zb].dmin,0,0) = nodata;
      val(_terrain[zb].dmax,0,0) = nodata;
    }
  }
}

void terrain (scalar zb, const char * name)
{
  if (zb >= _nterrain) {
    _terrain = realloc (_terrain, sizeof (Terrain)*(zb + 1));
    _nterrain = zb + 1;
  }
  _terrain[zb].kdt = kdt_new();
  if (kdt_open (_terrain[zb].kdt, name)) {
    fprintf (stderr, "terrain: could not open terrain database '%s'\n", name);
    exit (1);
  }
  zb.refine = refine_terrain;
  scalar n = new scalar;
  n.refine = refine_none;
  _terrain[zb].n = n;
  scalar dmin = new scalar;
  dmin.refine = refine_none;
  _terrain[zb].dmin = dmin;
  scalar dmax = new scalar;
  dmax.refine = refine_none;
  _terrain[zb].dmax = dmax;
}
