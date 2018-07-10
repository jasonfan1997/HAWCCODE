#!/bin/sh
indir=$1
outdir=$2
space=" "
prefix="./runone.sh"


for i in `find $1 -name "reco*" -type f`; do
    output=$prefix$space$i$space$2
	#$output
	echo $output >>temp.txt
done

