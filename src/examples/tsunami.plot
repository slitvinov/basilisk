! awk '{ if ($1 == "file:") file = $2; else print $0 > file; }' < tsunami.out

set term pngcairo enhanced size 700,700 font ",8"

set palette defined ( 0 0 0 0.5647, 0.125 0 0.05882 1, 0.25 0 0.5647 1, \
    	              0.375 0.05882 1 0.9333, 0.5 0.5647 1 0.4392, \
		      0.625 1 0.9333 0, 0.75 1 0.4392 0, \
		      0.875 0.9333 0 0, 1 0.498 0 0 )

unset key
set size ratio -1
set xrange [67:105]
set yrange [-10:25]
set pm3d map
set logscale cb
set cbrange [0.1:10]
splot 't-600' u 1:2:($3 > 1e-3 ? $5 : 1e1000)

! mogrify -trim +repage tsunami.png

reset
set term pngcairo enhanced size 625,800 font ",8"
set output 'tsunami_gauges.png'
set multiplot layout 5,1 scale 1,1.1
set xrange [3:8]
set key bottom right
set title 'Hanimaadhoo, Maldives'
plot 'hani' u ($1/60.):2 w l t 'modelled', \
     'hanires.txt' u 1:($2/100.) w lp t 'observed'
unset key
set title 'Male, Maldives'
plot 'male' u ($1/60.):2 w l t 'modelled', \
     'maleres.txt' u 1:($2/100.) w lp t 'observed'
set title 'Gan, Maldives'
plot 'gana' u ($1/60.):2 w l t 'modelled', \
     'ganares.txt' u 1:($2/100.) w lp t 'observed'
set title 'Diego Garcia'
plot 'dieg' u ($1/60.):2 w l t 'modelled', \
     'diegres.txt' u 1:($2/100.) w lp t 'observed'
set title 'Columbo, Sri Lanka'
set xrange [2.5:8]
plot 'colo' u ($1/60.):2 w l t 'modelled', \
     'colores.txt' u 1:($2/100.) w lp t 'observed'
unset multiplot

! rm -f t-*
