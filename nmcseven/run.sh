#!/bin/sh
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
pre="Aftercut_%s_"
ppre="Aftercut_"
suf="_.xcd"
und="_"
slash="/"

for i in $list
do

inputfilename="${i##*/}"
#echo $inputfilename
temp=${inputfilename/.xcd/en.xcd}

outfile=$pre$inputfilename
outname=${inputfilename/.xcd/add.xcd}
t="14"
test="/lustre/hawc04/scr01_backup/userspace/jasonfan/crab/test/14/"
if [ ! -e $test$ppre$t$und$inputfilename ]
then
	aerie-apps-add-energies --in-files $1$slash$inputfilename -o $2$temp
	python addvariable.py $2$temp $2
	./add.exe $2$temp $2${temp/en.xcd/add.xcd} $2${inputfilename/.xcd/en.csv}
	rm $2${inputfilename/.xcd/en.csv}
	rm $2$temp

	aerie-apps-scrappy-cut-dataset --cuts /data/disk01/home/jasonfan/nmcseven/cutmcnew --outfileFormat $2$outfile < $2$outname
	rm $2$outname
	#for i in `find $test -name "Aftercut_*_$inputfilename" -type f`; do
	for i in {0..119}
	do 
		if [ -e $2$ppre$i$und$inputfilename ]; then
			dir=$2$i$slash
			event=`xcdf count $2$ppre$i$und$inputfilename`
				if [ "$event" -eq "0" ]
				then
				rm $2$ppre$i$und$inputfilename
				else
					if [ -e $dir$ppre$i$und$inputfilename ]
					then
					
					rm $2$ppre$i$und$inputfilename
					else
					mv -n $2$ppre$i$und$inputfilename $dir
					
					
					#rm $2$ppre$i$und$inputfilename
					fi
				fi
			fi
	done
fi
done

for i in {0..119}
do
echo $i
lss=`find $2$i$slash -name $ppre$i$und"*"$run"*" |sort`
if [ "$lss" ]
then
xcdf select "1" -o $2$i$slash$run$und$i".xcd" $lss
rm $lss
fi

done
