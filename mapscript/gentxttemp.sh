indir=$1
outdir=$2
space=" "
prefix="makemaptemp.sh"
batch="sbatch --time=48:00:00 --mem-per-cpu=3000mb -n 1 -c 1 "
map="map/"
#rm temp.txt

for i in `find $1 -maxdepth 1 -type f -name "*_N1024.fits.gz"`; do
	input="${i##*/}"
	temp=${input/"_N1024.fits.gz"/}
	if [ ! -f $1$temp".fits.gz" ]
	then
	#rad="`head -n $i angle.csv | tail -1`"
    output=$batch$prefix$space$1$input$space$temp$space$1
	#$output
	echo $output >>submap.sh
	fi
done
