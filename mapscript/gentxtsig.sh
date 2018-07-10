indir=$1
outdir=$2
space=" "
prefix="makesigmap.sh"
batch="sbatch --time=48:00:00 --mem-per-cpu=4000mb -n 1 -c 1 "
map="map/"
#rm temp.txt
for i in `find $1 -maxdepth 1 -type f -name "*_N1024_smooth.fits.gz"`; do
	input="${i##*/}"
	temp=${input/"_N1024_smooth.fits.gz"/}
	#echo $temp
	if [ ! -f $1$temp".png" ]
	then
	#rad="`head -n $i angle.csv | tail -1`"
    output=$batch$prefix$space$1$input$space$temp$space$1$space$2
	#$output
	echo $output >>submap.sh
	#rm $1$temp".fits.gz"
	fi
	
done