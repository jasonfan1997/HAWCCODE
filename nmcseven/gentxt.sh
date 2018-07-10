#!/bin/sh
indir=$1
outdir=$2
space=" "
prefix="./run.sh"


for i in `find $1 -name "run*" -type d`; do
	number="`ls $i |wc -l`"
    output=$prefix$space$i$space$2$space$number
	#$output
	echo $output >>temp.txt
done

