indir=$1
outdir=$2
space=" "
prefix="makemap.sh"
batch="sbatch --time=48:00:00 --mem-per-cpu=4000mb -n 1 -c 1 "
map="map/"
#rm temp.txt
for i in `find $1 -type d`; do
	input="${i##*/}"
	if [ "$(ls -A $1$input)" ]
	then
	#rad="`head -n $i angle.csv | tail -1`"
    output=$batch$prefix$space$1$input$space$input$space$1$map
	#$output
	echo $output >>submap.sh
	fi
done

