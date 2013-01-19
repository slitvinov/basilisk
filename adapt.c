int coarsen_wavelet (Data * m, int n, var w, double max)
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
	  finecell(k,l).flags &= !leaf;
	  finecell(k,l).flags |= inactive;
	  /* update neighborhood */
	  for (int o = -GHOSTS; o <= GHOSTS; o++)
	    for (int p = -GHOSTS; p <= GHOSTS; p++)
	      finecell(k+o,l+p).neighbors--;
	}
      nc++;
    }
    /* propagate the error to coarser levels */
    coarse(w,0,0) = fabs(coarse(w,0,0)) + error;
  } end_foreach_fine_to_coarse ();
  return nc;
}

int flag_halo_cells (Data * m, int n)
{
  int nh = 0;
  /* from the bottom up */
  foreach_cell_post (m, n, cell.neighbors > 0 || !(cell.flags & inactive)) {
    if (cell.flags & inactive) {
      /* inactive and neighbors > 0 => this is a halo cell */
      cell.flags |= halo;
      /* propagate to parent */
      parentcell.flags |= halo;
      nh++;
    }
    else if (cell.flags & halo)
      /* propagate to parent */
      parentcell.flags |= halo;
  } end_foreach_cell_post();
  return nh;
}
