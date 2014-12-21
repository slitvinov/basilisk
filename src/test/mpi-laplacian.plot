reset
set multiplot layout 3,2

set ylabel 'Speed (points/sec/core)'
set grid

set yrange [*:*]
set logscale
set title "refine"
plot       '< grep refine res-8-8' u 1:($4/$1) w lp t '8 levels', \
	   '< grep refine res-8-9' u 1:($4/$1) w lp t '9 levels', \
	   '< grep refine res-8-13' u 1:($4/$1) w lp t '10 levels', \
	   '< grep refine res-8-11' u 1:($4/$1) w lp t '11 levels'

unset logscale
set logscale x
set yrange [0:]

do for [name in "cos laplacian sum restriction poisson"] {
   set title name
   plot '< grep '.name.' res-8-8' u 1:($4/$1) w lp t '8 levels', \
   	   '< grep '.name.' res-8-9' u 1:($4/$1) w lp t '9 levels', \
	   '< grep '.name.' res-8-13' u 1:($4/$1) w lp t '13 levels', \
	   '< grep '.name.' res-8-11' u 1:($4/$1) w lp t '11 levels'
}

unset multiplot

pause -1

reset
set multiplot layout 3,2

set ylabel 'Time (sec)'
set grid
set key top left
set yrange [0:]
set logscale x

do for [name in "refine cos laplacian sum restriction poisson"] {
   set title name
   plot '< grep '.name.' res-8-13' u 1:($3*$1) w lp t 'real', \
        '< grep '.name.' res-8-13' u 1:($6*$1) w lp t 'mpi (min)', \
        '< grep '.name.' res-8-13' u 1:($7*$1) w lp t 'mpi (avg)', \
        '< grep '.name.' res-8-13' u 1:($8*$1) w lp t 'mpi (max)', \
        '< grep '.name.' res-8-13' u 1:(($3-$7)*$1) w lp t 'real - mpi
}

unset multiplot
