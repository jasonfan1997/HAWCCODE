#!/bin/sh
indir=$1

space=" "
prefix="sbatch --time=48:00:00 --mem-per-cpu=400mb -n 1 -c 1 --output temp.out  all.sh "
for i in `find $1 -name "node*" -type f`; do
    output=$prefix$space$i
	#$output
	echo $output >>sub.sh
done
