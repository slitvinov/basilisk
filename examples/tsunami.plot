! awk '{ if ($1 == "file:") file = $2; else print $0 > file; }' < tsunami.out

set term pngcairo enhanced size 800,800 font ",11"

unset key
set size ratio -1
set xrange [67:105]
set yrange [-10:25]
set pm3d map
set logscale cb
set cbrange [0.1:10]
splot 't-600' u 1:2:($3 > 1e-3 ? $5 : 1e1000)

! rm -f t-?
