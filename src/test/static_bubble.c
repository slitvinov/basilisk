/**
# Soluble gas diffusing from a static bubble

This is the example discussed in section 3.3.1 of [Farsoiya et al.,
2020](#farsoiya2020).

<center>
<table>
<tr>
<td>![](static_bubble/final.png){ width="400px" }</td>
<td>![](static_bubble/p001cbt.png){ width="400px" }</td>
<td>![](static_bubble/p001clt.png){ width="400px" }</td>
</tr>
<tr>
<td><center>Diffusion from bubble</center></td> 
<td><center>Concentration field inside bubble</center></td> 
<td><center>Concentration field outside bubble</center></td> 
</tr>
</table>
</center>

The concentration at a point inside and outside the bubble is compared with the analytical solution provided in [Farsoiya et al.,
2020](#farsoiya2020).

~~~pythonplot

import numpy as np
import matplotlib.pyplot as plt
	
plt.figure()
t,Gb,Gl = np.loadtxt('../ref_static_bubble',delimiter=' ',unpack=True);
plt.plot(t/40, Gb,'k',label='Analytical')

ts,cb9,cl9 = np.loadtxt('out',delimiter=' ',unpack=True)
plt.plot(ts/40,cb9,'k--',label=r'$d_0/\Delta x \approx 102$');
# plt.ylim(0.97,1.007)

plt.legend();
plt.xlabel(r'$t\; \mathscr{D}_l/d_0^2$')
plt.ylabel(r'$c_b/c_{b0}$')
plt.tight_layout()

plt.savefig('p001cbt.png')

plt.figure()

plt.plot(t/40, Gl,'k',label='Analytical')

plt.plot(ts/40,cl9,'k--',label=r'$d_0/\Delta x \approx 102$');
# plt.ylim(0.97,1.007)

plt.legend();
plt.xlabel(r'$t\; \mathscr{D}_l/d_0^2$')
plt.ylabel(r'$c_l/c_{b0}$')
plt.tight_layout()

plt.savefig('p001clt.png')
plt.figure()

plt.savefig(' ')

~~~

## References

~~~bib
@article{farsoiya2020,
  title = {Bubble mediated gas transfer of diluted component in turbulence},
  author = {P. K. Farsoiya and S. Popinet and L. Deike},
  journal = {Journal of Fluid Mechanics},
  year = {2020},
  note = {submitted}
}
~~~
*/

#include "axi.h"
#include "navier-stokes/centered.h"
#include "two-phase.h"
#include "navier-stokes/conserving.h"
#include "tension.h"
#include "reduced.h"
#include "henry.h"
#include "view.h"

scalar c[], * stracers = {c};
double bubble_radius = 1.;
double box_size = 10.;
double conc_liq1 = 0, conc_gas1 = 1.;
	
int MAXLEVEL;

int main (int argc, char **argv)
{
  size (box_size);
	
  MAXLEVEL = 9;
  N = 1 << MAXLEVEL;


  rho1 = 1.;
  rho2 = 0.01;
  c.alpha = 0.001;
  TOLERANCE = 1e-4;
	
 
  c.D1 = 0.1;
  c.D2 = 1.;
    run();
  
}

event init (t = 0)
{

  fraction (f, - (sq(bubble_radius) - sq(x - box_size*0.5) - sq(y)));

  foreach()
    c[] = conc_liq1*f[] + conc_gas1*(1. - f[]);
}

	
event extract (t = 0; t += 1.; t <= 48.)
{		
  if (i == 0)
    fprintf (stdout, "#t ci co\n");
  fprintf(stdout,"%e %.12e %.12e\n",t,interpolate(c,5.0,0.5),interpolate(c,5.0,1.5));
}

event pictures (t = end)
{
  char name[80];
  dump(file="end");
  view (fov = 9, 
	tx = -0.5,
	width = 400, height = 400);  
  squares ("c", spread = 1, linear = true, map = cool_warm);
  draw_vof ("f");
  mirror ({0,1}) {
    squares ("c", spread = 1, linear = true, map = cool_warm);
    draw_vof ("f");
  }
  sprintf (name, "final.png");
  save (name);
}
