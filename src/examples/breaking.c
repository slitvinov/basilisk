/**
# 3D breaking Stokes wave (multilayer solver)

A steep, 3D, third-order Stokes wave is unstable and breaks. This is
the 3D equivalent of this [test case](/src/test/stokes.c). The
bathymetry is given by
$$
z_b(y) = - 0.5 + \sin(\pi y)/4
$$
i.e. it is shallower toward the back of the domain which causes the
wave to break earlier there.

![Animation of the free-surface. The surface is coloured according to
 the $x$-component of the surface velocity.](breaking/movie.mp4)(
 width=100% )

The solution is obtained using the layered model and demonstrates its
robustness and a degree of realism even for this complex case. Note
also the interesting longitudinal "scars". */

#include "grid/multigrid.h"
#include "layered/hydro.h"
#include "layered/nh.h"
#include "layered/remap.h"
#include "layered/perfs.h"

/**
The initial conditions are given by the wave steepness $ak$ and the
Reynolds number $Re=c\lambda/\nu$ with $c$ the phase speed of the
gravity wave and $\lambda$ its wavelength. */

double ak = 0.33;
double RE = 40000.;

#define k_  (2.*pi)
#define h_   0.5
#define g_   1.
#define T0  (k_/sqrt(g_*k_))

/**
The domain is periodic in $x$ and resolved using 256$^2$
points and 60 layers. */

int main()
{
  origin (-L0/2., -L0/2.);
  periodic (right);
  N = 256;
  nl = 60;
  G = g_;
  nu = 1./RE;
  run();
}

/**
The intial conditions for the free-surface and velocity are given by
the third-order Stokes solution. */

#include "../test/stokes.h"

event init (i = 0)
{
  foreach() {
    zb[] = - 0.5 + sin(pi*y)/4.;
    double H = wave(x, 0) - zb[];
    double z = zb[];
    vector u;
    scalar h, w;
    for (h,u,w in hl,ul,wl) {
      h[] = H/nl;
      z += h[]/2.;
      u.x[] = u_x(x, z);
      w[] = u_y(x, z);
      z += h[]/2.;
    }
  }
}

/**
We log the evolution of the kinetic and potential energies.

~~~gnuplot Evolution of the kinetic, potential and total energy
set xlabel 't/T0'
plot [0:6]'log' u 1:2 w l t 'kinetic', '' u 1:($3+0.0007) w l t 'potential', \
     '' u 1:(($2+$3+0.0007)/2.) w l t 'total/2'
~~~
*/

event logfile (i++)
{
  double ke = 0., gpe = 0.;
  foreach (reduction(+:ke) reduction(+:gpe)) {
    scalar h, w;
    vector u;
    double zc = zb[];
    for (h,w,u in hl,wl,ul) {
      double norm2 = sq(w[]);
      foreach_dimension()
	norm2 += sq(u.x[]);
      ke += norm2*h[]*dv();
      gpe += (zc + h[]/2.)*h[]*dv();
      zc += h[];
    }
  }
  fprintf (ferr, "%g %g %g\n", t/T0, ke/2., g_*gpe + 0.14);
}

/**
Note that the movie generation below is very expensive. */

#if 1
event movie (i += 5; t <= 6.*T0)
{
  static FILE * fp = popen ("gnuplot", "w");
  if (i == 0)
    fprintf (fp, "set term pngcairo font ',9' size 1024,768;"
	     "unset key\n"
	     "set pm3d interpolate 4,4 lighting specular 0.6\n"
	     "set zrange [-0.35:0.15]\n"
	     "set cbrange [-0.15:0.6]\n"
	     "set xlabel 'x'\n"
	     "set ylabel 'y'\n"
	     );
  fprintf (fp,
	   "set output 'plot%04d.png'\n"
	   "set title 't = %.2f T0'\n"
	   "splot '-' u 1:2:3:4 w pm3d\n",
	   i/3, t/T0);
  scalar H[];
  foreach() {
    H[] = zb[];
    for (scalar h in hl)
      H[] += h[];
  }
  boundary ({H});
  scalar ux = ul[nl-1].x;
  output_field ({H,ux}, fp, linear = true);
  fprintf (fp, "e\n\n");
  fflush (fp);
}

event moviemaker (t = end)
{
  system ("for f in plot*.png; do convert $f ppm:-; done | ppm2mp4 movie.mp4");
}
#endif

/**
## Parallel run

The simulation was run in parallel on the [Jean Zay
machine](http://www.idris.fr/jean-zay/) on 64 cores, using this script

~~~bash
local% qcc -source -D_MPI=1 breaking.c
local% scp _breaking.c user@jean-zay.idris.fr:
~~~

~~~bash
#!/bin/bash
#SBATCH --job-name=breaking       # nom du job
#SBATCH --ntasks=64              # Nombre total de processus MPI
#SBATCH --ntasks-per-node=32       # Nombre de processus MPI par noeud
# /!\ Attention, la ligne suivante est trompeuse mais dans le vocabulaire
# de Slurm "multithread" fait bien référence à l'hyperthreading.
#SBATCH --hint=nomultithread       # 1 processus MPI par coeur physique (pas d'hyperthreading)
#SBATCH --time=2:00:00            # Temps d’exécution maximum demande (HH:MM:SS)
#SBATCH --output=breaking%j.out  # Nom du fichier de sortie
#SBATCH --error=breaking%j.out   # Nom du fichier d'erreur (ici commun avec la sortie)
#SBATCH --mail-type=BEGIN,END,FAIL
#SBATCH --mail-user=popinet@basilisk.fr
 
# on se place dans le répertoire de soumission
cd ${SLURM_SUBMIT_DIR}
 
# nettoyage des modules charges en interactif et herites par defaut
module purge
 
# chargement des modules
module load intel-all/19.0.4 gnuplot
 
# echo des commandes lancées
set -x

# compilation
mpiicc -Wall -O2 -std=c99 -I$HOME -I$HOME/local/include -L$HOME/local/lib _breaking.c -o breaking -lppr -lgfortran -lm
 
# exécution du code
export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:$HOME/local/lib
rm -f plot*.png
srun ./breaking 2> log > out
~~~

The number of timesteps was 4334 and the runtime was 73 minutes with
movie generation and 13 minutes without, corresponding to a
computational speed of 338 000 point.timestep/sec/core (on 64 cores). */
