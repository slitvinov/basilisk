! awk '{ if ($1 == "file:") file = $2; else print $0 > file; }' < conical.out
 
set term pngcairo enhanced size 1600,400 font ",8"

unset key
set hidden3d front
unset xtics
unset ytics
set xyplane at 0
unset colorbox

dry=1e-4
set view 28,32

set multiplot layout 1,3 scale 1.2,1.2
splot './conical-0' u 1:2:($3>dry?$3+$4:1e1000) w pm3d, \
      './conical-0' u 1:2:4 every 2:2 w l, \
      './conical-0' u 1:2:($3>dry?$3+$4:1e1000) w l lt -2
splot './conical-1' u 1:2:($3>dry?$3+$4:1e1000) w pm3d, \
      './conical-0' u 1:2:4 every 2:2 w l, \
      './conical-1' u 1:2:($3>dry?$3+$4:1e1000) w l lt -2
splot './conical-2' u 1:2:($3>dry?$3+$4:1e1000) w pm3d, \
      './conical-0' u 1:2:4 every 2:2 w l, \
      './conical-2' u 1:2:($3>dry?$3+$4:1e1000) w l lt -2
unset multiplot

reset
set term pngcairo enhanced size 640,640
set output 'conical_runup.png'
set size ratio -1
unset surface
unset key
set contour base
set cntrparam levels discrete 1e-2
set table 'runup.dat'
splot './conical-3' u 1:2:5 w l
unset table
set grid
plot [10:16][11:17]'./runup.dat' w l

! rm -f conical-?
