int coarsen_wavelet (void * m, int n, var w, double max)
{
  int nc = 0;
  foreach_fine_to_coarse (m, n) {
    double error = 0.;
    foreach_fine(k,l) {
      double e = fabs(fine(w,k,l));
      if (e > error)
	error = e;
    }
    if (error < max) {
      /* coarsen */
      cell.flags |= leaf;
      for (int k = 0; k < 2; k++)
	for (int l = 0; l < 2; l++) {
	  child(k,l).flags &= !leaf;
	  child(k,l).flags |= inactive;
#if 1
	  /* trash the data */
	  for (int v = 0; v < sizeof(Data)/sizeof(double); v++)
	    ((double *)&child(k,l).d)[v] = -1e30;
#endif
	  /* update neighborhood */
	  for (int o = -GHOSTS; o <= GHOSTS; o++)
	    for (int p = -GHOSTS; p <= GHOSTS; p++)
	      child(k+o,l+p).neighbors--;
	}
      nc++;
    }
    /* propagate the error to coarser levels */
    val(w,0,0) = fabs(val(w,0,0)) + error;
  } end_foreach_fine_to_coarse ();
  return nc;
}

int flag_halo_cells (void * m, int n)
{
  int nh = 0;
  /* from the bottom up */
  foreach_cell_post (m, n, cell.neighbors > 0 || !(cell.flags & inactive)) {
    if (cell.flags & inactive) {
      /* inactive and neighbors > 0 => this is a halo cell */
      cell.flags |= halo;
      /* propagate to parent */
      parent.flags |= halo;
      nh++;
    }
    else if (cell.flags & halo)
      /* propagate to parent */
      parent.flags |= halo;
  } end_foreach_cell_post();
  return nh;
}

void update_halos (void * m, int n, var start, var end)
{
  /* breadth-first traversal of halos from coarse to fine */
  for (int l = 1; l <= mlevel(n); l++)
    /* fixme: need to stop before mlevel(n) if there an't any halos */
    foreach_cell(m,n) {
      if (!(cell.flags & halo))
	continue; /* no more halos, skip the rest of this branch */
      if (level == l) {
	if (cell.flags & inactive)
	  /* bilinear interpolation from coarser level */
	  for (var a = start; a <= end; a += sizeof(double))
	    val(a,0,0) = 
	      (9.*coarse(a,0,0) + 
	       3.*(coarse(a,childx,0) + coarse(a,0,childy)) + 
	       coarse(a,childx,childy))/16.;
	continue; /* already at level l, skip the deeper branches */
      }
    } end_foreach_cell();
}
