# checks the consistency of send/receive buffers

npe=3 # the number of procs - 1
for i in `seq 0 1 $npe`; do
    for j in `seq 0 1 $npe`; do
	if test -f rcv-$i-$j -a -f snd-$j-$i; then
	    echo rcv-$i-$j snd-$j-$i
	    diff rcv-$i-$j snd-$j-$i
	fi
    done
done
