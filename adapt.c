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
#if 1
      /* trash the data just to make sure it's never touched */
      for (int v = 0; v < sizeof(Data)/sizeof(double); v++)
	((double *)&child(k,l).d)[v] = -1e30;
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
    val(w,0,0) = fabs(val(w,0,0)) + error;
  }
  return nc;
}

/* fixme: this needs to be merged with refine() in quadtree.c */
Point refine_cell (Point point, var start, var end)
{
  QUADTREE_VARIABLES;
  alloc_children();
  cell.flags &= ~leaf;
  for (int k = 0; k < 2; k++)
    for (int l = 0; l < 2; l++) {
      assert(!(child(k,l).flags & active));
      child(k,l).flags |= (active | leaf);
      /* update neighborhood */
      for (int o = -GHOSTS; o <= GHOSTS; o++)
	for (int p = -GHOSTS; p <= GHOSTS; p++)
	  child(k+o,l+p).neighbors++;
      /* bilinear interpolation from coarser level */
      for (var a = start; a <= end; a += sizeof(double)) /* for each variable */
	fine(a,k,l) = 
	  (9.*val(a,0,0) + 3.*(val(a,2*k-1,0) + val(a,0,2*l-1)) + val(a,2*k-1,2*l-1))/16.;
    }
  return point;
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
    if (fabs(val(w,0,0)) >= max) {
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

void update_halos (Quadtree * quadtree, int depth, var start, var end)
{
  if (depth < 0)
    depth = quadtree->depth; /* fixme */
  /* breadth-first traversal of halos from coarse to fine */
  for (int l = 0; l <= depth; l++)
    foreach_cell (quadtree) {
      if (!(cell.flags & halo))
      	continue; /* no more halos, skip the rest of this branch */
      if (level == l) {
	if (!(cell.flags & active))
	  /* bilinear interpolation from coarser level */
	  for (var a = start; a <= end; a += sizeof(double)) /* for each variable */
	    val(a,0,0) = 
	      (9.*coarse(a,0,0) + 
	       3.*(coarse(a,childx,0) + coarse(a,0,childy)) + 
	       coarse(a,childx,childy))/16.;
	continue; /* already at level l, skip the deeper branches */
      }
    }
}

#define foreach_halo(quadtree) foreach_cell(quadtree) { \
  if (!(cell.flags & halo))		      \
    continue;				      \
  else if (!(cell.flags & active)) {
#define end_foreach_halo()  }} end_foreach_cell();
