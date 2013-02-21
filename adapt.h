void coarsen_cell (Point point)
{
  QUADTREE_VARIABLES;
  /* coarsen */
  free_children();
  cell.flags |= leaf;
  for (int k = 0; k < 2; k++)
    for (int l = 0; l < 2; l++) {
      child(k,l).flags &= ~leaf;
      child(k,l).flags &= ~active;
#if TRASH
      /* trash the data just to make sure it's never touched */
      for (var v = 0; v < nvar; v++)
	fine(v,k,l) = undefined;
#endif
      /* update neighborhood */
      for (int o = -GHOSTS; o <= GHOSTS; o++)
	for (int p = -GHOSTS; p <= GHOSTS; p++)
	  child(k+o,l+p).neighbors--;
    }
}

int coarsen_function (Quadtree * quadtree, int (* func) (Point p))
{
  int nc = 0;
  foreach_fine_to_coarse (quadtree)
    if ((*func) (point)) {
      coarsen_cell (point);
      nc++;
    }
  return nc;
}

int coarsen_wavelet (Quadtree * quadtree, var w, double max)
{
  int nc = 0;
  foreach_fine_to_coarse (quadtree) {
    double error = 0.;
    for (int k = 0; k < 2; k++)
      for (int l = 0; l < 2; l++) {
	double e = fabs(fine(w,k,l));
	if (e > error)
	  error = e;
      }
    if (error < max) {
      coarsen_cell(point);
      nc++;
    }
    /* propagate the error to coarser levels */
    w(0,0) = fabs(w(0,0)) + error;
  }
  return nc;
}

int refine_function (Quadtree * quadtree, var start, var end,
		     int (* func) (Point p, void * data), 
		     void * data)
{
  int nf = 0;
  foreach_leaf (quadtree)
    if ((*func) (point, data)) {
      point = refine_cell (point, start, end);
      nf++;
    }
  return nf;
}

int refine_wavelet (Quadtree * quadtree, var start, var end,
		    var w, double max)
{
  int nf = 0;
  foreach_leaf (quadtree)
    if (fabs(w(0,0)) >= max) {
      point = refine_cell (point, start, end);
      nf++;
    }
  return nf;
}

int flag_halo_cells (Quadtree * quadtree)
{
  int nh = 0;

  /* reset old halos first */
  foreach_cell (quadtree) {
    if (!(cell.flags & halo))
      continue;
    else 
      cell.flags &= ~halo;
  }

  /* from the bottom up */
  foreach_cell_post (quadtree, cell.neighbors > 0 || (cell.flags & active)) {
    if (!(cell.flags & active)) {
      /* inactive and neighbors > 0 => this is a halo cell */
      cell.flags |= halo;
      /* propagate to parent */
      parent.flags |= halo;
      nh++;
    }
    else if ((cell.flags & halo) && level > 0)
      /* propagate to parent */
      parent.flags |= halo;
  }

  return nh;
}

#define foreach_halo(quadtree) foreach_cell(quadtree) { \
  if (!(cell.flags & halo))		      \
    continue;				      \
  else if (!(cell.flags & active)) {
#define end_foreach_halo()  }} end_foreach_cell();

/* breadth-first traversal of halos from coarse to fine */
#define foreach_halo_coarse_fine(quadtree,depth1)    {			\
  int _depth = depth1 < 0 ? depth(quadtree) : depth1;			\
  for (int _l = 0; _l <= _depth; _l++)                                  \
    foreach_cell (quadtree) {                                           \
      if (!(cell.flags & halo))                                         \
      	continue; /* no more halos, skip the rest of this branch */     \
      if (level == _l) {                                                \
	if (!(cell.flags & active))
#define end_foreach_halo_coarse_fine()                                  \
	continue; /* already at level l, skip the deeper branches */    \
      }                                                                 \
    } end_foreach_cell();				                \
}

void update_halo (Quadtree * quadtree, int depth, var start, var end)
{
  foreach_halo_coarse_fine (quadtree, depth)
    /* bilinear interpolation from coarser level */
    for (var v = start; v <= end; v++)
      val(v,0,0) = (9.*coarse(v,0,0) + 
		    3.*(coarse(v,childx,0) + coarse(v,0,childy)) + 
		    coarse(v,childx,childy))/16.;
}

void update_halo_u_v (Quadtree * quadtree, int depth, var u, var v)
{
  foreach_halo_coarse_fine (quadtree, depth) {
    /* linear interpolation from coarser level */
    u[] = coarse(u,0,0) + childy*(coarse(u,0,1) - coarse(u,0,-1))/8.;
    v[] = coarse(v,0,0) + childx*(coarse(v,1,0) - coarse(v,-1,0))/8.;
  }
}