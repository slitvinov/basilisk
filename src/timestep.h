double timestep (const vector u)
{
  double dtmax = DT/CFL;
  foreach(reduction(min:dtmax))
    foreach_dimension()
      if (u.x[] != 0.) {
	double dt = delta/fabs(u.x[]);
	if (dt < dtmax) dtmax = dt;
      }
  return dtmax*CFL;
}
