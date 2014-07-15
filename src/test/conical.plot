! awk '{ if ($1 == "file:") file = $2; else print $0 > file; }' < out
 
set term @PNG enhanced size 900,1600 font ",8"

unset key
unset xtics
unset ytics

dry=1e-4

set size ratio -1
set pm3d
set pm3d map
# set contour base
set palette defined ( 0 0 0 0.5647, 0.125 0 0.05882 1, 0.25 0 0.5647 1,\
     0.375 0.05882 1 0.9333, 0.5 0.5647 1 0.4392, 0.625 1 0.9333 0, \
     0.75 1 0.4392 0, 0.875 0.9333 0 0, 1 0.498 0 0 )

set multiplot layout 4,2 scale 1.35,1.35
splot './conical-0' u 1:2:($3>dry?$3+$4:1e1000)
splot './level-0'
splot './conical-1' u 1:2:($3>dry?$3+$4:1e1000)
splot './level-1'
splot './conical-2' u 1:2:($3>dry?$3+$4:1e1000)
splot './level-2'
splot './conical-3' u 1:2:($3>dry?$3+$4:1e1000)
splot './level-3'
unset multiplot

reset
set term @PNG enhanced size 640,640
set output 'runup.png'
set size ratio -1
unset surface
unset key
set contour base
set cntrparam levels discrete 1e-2,1e-3
set table 'runup.dat'
splot './conical-4' u 1:2:5 w l
unset table
set grid
plot [10:16][11:17]'./runup.dat' w l

reset
set term @PNG enhanced size 600,800
set output 'gauges.png'
set multiplot layout 5,1 scale 1.,1.
set xrange [3:20]
t0=22.
h0=0.32
plot '../ts2b.txt' u ($1-t0):($4-0.001) pt 7 ps 0.25 t '', \
     '../conicalsv/WG3' u 1:($2-h0) w l lt 4 t '', \
     './WG3' u 1:($2-h0) w l lw 2 lt 3 t 'WG3'
plot '../ts2b.txt' u ($1-t0):($6-0.0004) pt 7 ps 0.25 t '', \
     '../conicalsv/WG6' u 1:($2-h0) w l lt 4 t '', \
     './WG6' u 1:($2-h0) w l lw 2 lt 3 t 'WG6'
plot '../ts2b.txt' u ($1-t0):($7) pt 7 ps 0.25 t '', \
     '../conicalsv/WG9' u 1:($2-h0) w l lt 4 t '', \
     './WG9' u 1:($2-h0) w l lw 2 lt 3 t 'WG9'
plot '../ts2b.txt' u ($1-t0):($8-0.0017) pt 7 ps 0.25 t '', \
     '../conicalsv/WG16' u 1:($2-h0) w l lt 4 t '', \
     './WG16' u 1:($2-h0) w l lw 2 lt 3 t 'WG16'
plot '../ts2b.txt' u ($1-t0):($9+0.0015) pt 7 ps 0.25 t '', \
     '../conicalsv/WG22' u 1:($2-h0) w l lt 4 t '', \
     './WG22' u 1:($2-h0) w l lw 2 lt 3 t 'WG22'
unset multiplot

! rm -f conical-?
