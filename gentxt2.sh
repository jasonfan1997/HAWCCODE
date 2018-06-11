#export PATH="~/miniconda2/bin:$PATH";
#./run.sh /lustre/hawc01/hawcroot/data/hawc/data/subsets/crab-strip/2017/08/run007094/reco-crab-strip_run007094_00189.xcd ./data/
indir=$1
outdir=$2
space=" "
prefix="./run2.sh"
for i in `find $1 -name "reco*" -type f`; do
    output=$prefix$space$i$space$2
	#$output
	echo $output >>temp.txt
done