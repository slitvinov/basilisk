# checks the consistency of send/receive buffers

npe=4 # the number of procs - 1
for name in restriction prolongation halo-restriction; do
    for i in `seq 0 1 $npe`; do
	for j in `seq 0 1 $npe`; do
	    if test -f $name-rcv-$i-$j -a -f $name-snd-$j-$i; then
		echo $name-rcv-$i-$j $name-snd-$j-$i
		diff $name-rcv-$i-$j $name-snd-$j-$i
	    fi
	done
	# cat $name-rcv-$i-* | awk '{print $1,$2}' | sort > /tmp/$name-rcv-$i
	# awk '{print $1,$2}' mpi-$name-$i | sort > /tmp/mpi-$name-$i
	# echo $name-rcv-$i mpi-$name-$i
	# diff /tmp/$name-rcv-$i /tmp/mpi-$name-$i
    done
done
