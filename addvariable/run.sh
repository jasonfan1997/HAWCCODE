#!/bin/sh
#export PATH="~/miniconda2/bin:$PATH";
export LD_LIBRARY_PATH="~/miniconda2/lib:$LD_LIBRARY_PATH"
export PATH="~/miniconda2/envs/myconda/bin:$PATH"
eval `$HOME/utils/ape/current/ape sh externals`
eval `~/aerie/install/bin/hawc-config --env-sh`
export LD_LIBRARY_PATH=${CONDA_PREFIX}/lib:/data/disk01/software/hawc/ApeInstalled/x86_64/External/gsl/1.16/lib/:${LD_LIBRARY_PATH}
export CONFIG_HAWC=/data/disk01/home/jasonfan/config-hawc


input=$1
outdir=$2
pattern=$1"/reco*"
list=`ls $pattern`
run="${1##*/}"
pre="_%s"
ppre="Aftercut_"
suf=".xcd"
und="_"
slash="/"


input="${1##*/}"
inputrun=$input$".xcd"
outfile=$input$pre$suf
#echo $inputfilename

temp=${inputrun/.xcd/en.xcd}


outname=${inputrun/.xcd/add.xcd}

aerie-apps-add-energies --in-files $list -o $QTMPDIR$slash$temp
./addproba.exe $QTMPDIR$slash$temp $QTMPDIR$slash$outname
rm $QTMPDIR$slash$temp


aerie-apps-scrappy-cut-dataset --cuts /data/disk01/home/jasonfan/newmc/newmcmc --outfileFormat $2$outfile < $QTMPDIR$slash$outname
echo "finish cut"
rm $QTMPDIR$slash$outname

#for i in `find $test -name "Aftercut_*_$inputfilename" -type f`; do
for i in {0..119}
do 
	if [ -e $2$input$und$i$suf ]; then
		dir=$2$i$slash
		event=`xcdf count $2$input$und$i$suf`
			if [ "$event" -eq "0" ]
			then
			rm $2$input$und$i$suf
			else
				if [ -e $dir$input$und$i$suf ]
				then
				
				rm $2$input$und$i$suf
				else
				mv -n $2$input$und$i$suf $dir
				
				
				fi
			fi
		fi
done

