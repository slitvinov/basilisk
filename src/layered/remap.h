#include "ppr/ppr.h"

// int edge_meth = p1e_method, cell_meth = plm_method, cell_lim = null_limit;
int edge_meth = p3e_method, cell_meth = ppm_method, cell_lim = null_limit;
// int edge_meth = p5e_method, cell_meth = pqm_method, cell_lim = null_limit;

trace
void vertical_remapping (scalar * hl, scalar ** tracers)
{
  int nvar = list_len(tracers[0]) - 1, ndof = 1, npos = nl + 1;
  assert (tracers[0][0].i == hl[0].i);
  foreach() {
    // remap to sigma coordinates
    double H = 0.;
    for (scalar h in hl)
      H += h[];

    if (H > dry) {
      double zpos[npos], znew[npos];
      double fdat[nvar*nl], fnew[nvar*nl];
      zpos[0] = znew[0] = 0.;
      int l = 0;
      for (scalar h in hl) {
	zpos[l+1] = zpos[l] + h[];
	int i = nvar*l;
	scalar * list = tracers[l] + 1;
	for (scalar s in list)
	  fdat[i++] = s[];
#if 1
	h[] = beta[l]*H;
#else
	double z0 = - 0.1;
	if (l < nl/2)
	  h[] = (z0 - zb[])/(nl/2);
	else
	  h[] = (H + zb[] - z0)/(nl/2);
#endif
	znew[l+1] = znew[l] + h[];
	l++;
      }

#if 1      
      double H1 = 0.;
      for (scalar h in hl)
	H1 += h[];
      assert (fabs(H - H1) < 1e-8);
#endif
      
      my_remap (&npos, &npos, &nvar, &ndof, zpos, znew, fdat, fnew,
		&edge_meth, &cell_meth, &cell_lim);

      l = 0;
      for (scalar h in hl) {
	int i = nvar*l;
	scalar * list = tracers[l] + 1;
	for (scalar s in list)
	  s[] = fnew[i++];
	l++;
      }
    }
  }
  
  scalar * list = NULL;
  for (int l = 0; l < nl; l++) {
    scalar * ll = tracers[l];
    for (scalar s in ll)
      list = list_append (list, s);
  }
  boundary (list);
  free (list);
}

event viscous_term (i++) {
  if (nl > 1)
    vertical_remapping (hl, tracers);
}
