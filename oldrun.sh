#export PATH="~/miniconda2/bin:$PATH";
export LD_LIBRARY_PATH="~/miniconda2/lib:$LD_LIBRARY_PATH"
export PATH="~/miniconda2/envs/myconda/bin:$PATH"
eval `$HOME/utils/ape/current/ape sh externals`
eval `~/aerie/install/bin/hawc-config --env-sh`
export LD_LIBRARY_PATH=${CONDA_PREFIX}/lib:/data/disk01/software/hawc/ApeInstalled/x86_64/External/gsl/1.16/lib/:${LD_LIBRARY_PATH}
export CONFIG_HAWC=/data/disk01/home/jasonfan/config-hawc

input=$1
outdir=$2
inputfilename="${1##*/}"
python addvariable.py $1 $2
./add.exe $1 $2${inputfilename/.xcd/add.xcd} $2${inputfilename/.xcd/.csv}
rm $2${inputfilename/.xcd/.csv}
pre="Afterocut_%s_"
outfile=$pre$inputfilename
outname=${inputfilename/.xcd/add.xcd}
aerie-apps-scrappy-cut-dataset --cuts /data/disk01/home/jasonfan/crab/new --outfileFormat $2$outfile < $2$outname
#rm $2$outname
