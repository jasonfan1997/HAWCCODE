indir=$1
outdir=$2
space=" "
prefix="./run.sh"

rm temp.txt
for i in `find $1 -name "reco*" -type f`; do
	filename="${i##*/}"
	temp=${i#"$1"}
	test=${temp//$filename}
	final=$2$test
    output=$prefix$space$i$space$final
	#$output
	echo $output >>temp.txt
done

