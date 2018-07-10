#!/bin/sh
#export PATH="~/miniconda2/bin:$PATH";
export LD_LIBRARY_PATH="~/miniconda2/lib:$LD_LIBRARY_PATH"
export PATH="~/miniconda2/envs/myconda/bin:$PATH"
eval `$HOME/utils/ape/current/ape sh externals`
eval `~/aerie/install/bin/hawc-config --env-sh`
export LD_LIBRARY_PATH=${CONDA_PREFIX}/lib:/data/disk01/software/hawc/ApeInstalled/x86_64/External/gsl/1.16/lib/:${LD_LIBRARY_PATH}
export CONFIG_HAWC=/data/disk01/home/jasonfan/config-hawc
cd $1
afk="Aft*"
aerie-apps-make-hawc-maps --useJ2000 --dtMin_hr 0 --nSide 1024 -z $CONFIG_HAWC/reconstruction/crab-align/zenith-SPofficial.xml --input $afk -n $2 -p $3

