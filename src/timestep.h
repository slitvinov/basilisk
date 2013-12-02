double timestep (const face vector u)
{
  double dtmax = DT/CFL;
  foreach_face(reduction(min:dtmax))
    if (u.x[] != 0.) {
      double dt = Delta/fabs(u.x[]);
      if (dt < dtmax) dtmax = dt;
    }
  return dtmax*CFL;
}
