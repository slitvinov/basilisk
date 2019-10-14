double filter = 0.;

event viscous_term (i++)
{
  if (filter > 0.) {
    foreach()
      foreach_dimension() {
        double Hm = 0., H = 0., Hp = 0.;
	for (scalar h in hl)
	  Hm += h[-1], H += h[], Hp += h[1];
        if (Hm > dry && H > dry && Hp > dry) {
	  double Dp = eta[1] - eta[], Dm = eta[] - eta[-1];
	  if (Dp*Dm < 0. && ((eta[2] + eta[] - 2.*eta[1])*
			     (eta[-1] + eta[1] - 2.*eta[]) < 0. ||
			     (eta[-1] + eta[1] - 2.*eta[])*
			     (eta[-2] + eta[] - 2.*eta[-1]) < 0.)) {
	    double dp, dm;
	    if (fabs(Dp) > fabs(Dm)) {
	      dp = fabs(Dp);
	      dm = fabs(Dm);
	    }
	    else {
	      dp = fabs(Dm);
	      dm = fabs(Dp);
	    }
	    double d = min(dm, dp/2.);
	    double a = Dp > 0. ? 1. : -1.;
	    eta[] += min(dt/filter, 1.)*a*d;
	    double Hnew = eta[] - zb[];
	    if (Hnew > dry) {
	      for (scalar h in hl)
		h[] *= Hnew/H;
	    }
	    else {
	      for (int l = 0; l < nl; l++)
		for (scalar s in tracers[l])
		  s[] = 0.;
	    }
	  }
	}
      }
    boundary ({eta});
    boundary (hl);
    for (int l = 0; l < nl; l++)
      boundary (tracers[l]);
  }
}
