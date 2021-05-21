/**
# Horizontal viscosity

This is a naive implementation of horizontal viscosity for the
[multilayer solver](README). Metric terms linked to the slope of
layers are most probably missing (but horizontal metric terms are
taken into account). This also assumes that the divergence of the
horizontal velocity field in each layer is zero, which is obviously
not correct: non-diagonal viscous stresses would need to be added to
relax this assumption.

This approximates
$$
\begin{aligned}
\partial_t u & = \frac{\nu_H}{h} \nabla \cdot (h \nabla u) \\
\partial_t v & = \frac{\nu_H}{h} \nabla \cdot (h \nabla v)
\end{aligned}
$$
*/

double nu_H = 0.;

event viscous_term (i++)
{
  if (nu_H > 0.) {
    vector d2u[];
    foreach_layer() {
      foreach()
	foreach_dimension() {
	  scalar s = u.x;
	  double a = 0.;
	  foreach_dimension()
	    a += (hf.x[]*fm.x[]/(cm[-1] + cm[])*(s[-1] - s[]) +
		  hf.x[1]*fm.x[1]/(cm[1] + cm[])*(s[1] - s[]));
	  d2u.x[] = 2.*a/(cm[]*sq(Delta));
        }
      foreach()
	if (h[] > dry)
	  foreach_dimension()
	    u.x[] += dt*nu_H*d2u.x[]/h[];
    }
    boundary ((scalar *){u});
  }
}
