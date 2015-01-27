# generate results for Curie
cd '../curie'

# generate weak scaling curves
! bash weak.sh > weak

# change default color for line style 6 from "yellow" to "sea-green"
set style line 6 lw 1 lc rgb "sea-green"
set style increment user

set logscale
set grid
set xrange [1:16384]
set xtics 2

# Model for memory usage
cst(level)=2**(2*level)*12./1024**3
data(level)=2**(2*level)*8*20./1024**3
tot(level,np)=cst(level)/sqrt(np)+data(level)/np+np*7e-6+0.03

set xlabel "# of cores"
set ylabel "Memory/core (GB)"
set output 'memory.png'
plot [][0.01:]\
     for [i=10:15] '< sh table.sh poisson '.i u 1:($2/$1) t ''.i.' levels', \
     for [i=10:15] tot(i,x) t '' lt 1

set ylabel 'Time (sec)'

# model of computation time
tc(level,np)=2**(2.*level)/(np*1.5e6)
tcom(level,np)=2**(1.6*level)/(1e7*sqrt(np))+np**0.4/4e3
tt(level,np)=tc(level,np)+tcom(level,np)

set output 'poisson.png'
plot [][:100] \
     for [i=10:15] '< sh time.sh poisson '.i u 1:2 t ''.i.' levels', \
     for [i=10:15] tt(i,x) t '' lt 1, \
     'weak' u 1:2 w l t 'weak scaling'
     
set output 'poisson-mpi.png'
plot [][1e-3:]\
     for [i=10:15] '< sh time.sh poisson '.i u 1:3 t ''.i.' levels', \
     for [i=10:15] tcom(i,x) t '' lt 1

set output 'laplacian.png'
tc(level,np)=2**(2.*level)/(np*25e6)
tcom(level,np)=2**(1.6*level)/(1.5e8*sqrt(np))
plot [][1e-4:]\
     for [i=10:15] '< sh time.sh laplacian '.i u 1:2 t ''.i.' levels', \
     for [i=10:15] tt(i,x) t '' lt 1

set output 'laplacian-mpi.png'
plot [][1e-5:]\
     for [i=10:15] '< sh time.sh laplacian '.i u 1:3 t ''.i.' levels', \
     for [i=10:15] tcom(i,x) t '' lt 1

set output 'restriction.png'
tc(level,np)=2**(2.*level)/(np*30e6)
tcom(level,np)=2**(1.5*level)/(3.5e7*sqrt(np))+np**0.3/1.5e4
plot [][:10]\
     for [i=10:15] '< sh time.sh restriction '.i u 1:2 t ''.i.' levels', \
     for [i=10:15] tt(i,x) t '' lt 1

set output 'restriction-mpi.png'
plot [][:0.1]\
     for [i=10:15] '< sh time.sh restriction '.i u 1:3 t ''.i.' levels', \
     for [i=10:15] tcom(i,x) t '' lt 1
