indir=$1
outdir=$2
space=" "
prefix="./run.sh"

rm temp.txt
:'
for i in `find $1 -name "reco*" -type f`; do
	filename="${i##*/}"
	temp=${i#"$1"}
	test=${temp//$filename}
	final=$2$test
    output=$prefix$space$i$space$final
	#$output
	echo $output >>temp.txt
done
'
for i in `find $1 -name "reco*" -type f`; do
    output=$prefix$space$i$space$2
	#$output
	echo $output >>temp.txt
done
for i in {0..119}
do 
	dir=$2$i
	mkdir -p $dir
done
open_sem(){
    mkfifo pipe-$$
    exec 3<>pipe-$$
    rm pipe-$$
    local i=$1
    for((;i>0;i--)); do
        printf %s 000 >&3
    done
}
run_with_lock(){
    local x
    read -u 3 -n 3 x && ((0==x)) || exit $x
    (
    "$@" 
    printf '%.3d' $? >&3
    )&
}

N=8
open_sem $N
cat temp.txt | while read p
do
    run_with_lock $p
done 