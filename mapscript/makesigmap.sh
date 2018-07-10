#!/bin/sh
#export PATH="~/miniconda2/bin:$PATH";
export LD_LIBRARY_PATH="~/miniconda2/lib:$LD_LIBRARY_PATH"
export PATH="~/miniconda2/envs/myconda/bin:$PATH"
eval `$HOME/utils/ape/current/ape sh externals`
eval `~/aerie/install/bin/hawc-config --env-sh`
export LD_LIBRARY_PATH=${CONDA_PREFIX}/lib:/data/disk01/software/hawc/ApeInstalled/x86_64/External/gsl/1.16/lib/:${LD_LIBRARY_PATH}
export CONFIG_HAWC=/data/disk01/home/jasonfan/config-hawc
#cd $1

fits=".fits.gz"
empty=""
pre=$2
outname=$pre$fits
png=".png"
afk="Aft*"
#aerie-apps-make-hawc-maps --useJ2000 --dtMin_hr 0 --nSide 1024 -z $CONFIG_HAWC/reconstruction/crab-align/zenith-SPofficial.xml --input $1 --itype crawl -n $2 -p $3

#aerie-apps-HealpixSigFluxMap -i $3$pre"_N1024.fits.gz" -b $pre -d /lustre/hawc04/scr01_backup/userspace/jasonfan/crab/response -o $3$outname --negFlux --window 83.6332 22.0145 4 4 -s SimplePowerLaw,3.5e-11,2.63
#aerie-apps-recalculate-bkg -i $3$pre"_N1024.fits.gz" -o $3$pre"_N1024_smooth.fits.gz" --smooth 0.05,0
aerie-apps-HealpixSigFluxMap -i $1 -b $pre -d /lustre/hawc04/scr01_backup/userspace/jasonfan/crab/response -o $3$outname --negFlux --window 83.6332 22.0145 4 4 -s SimplePowerLaw,3.5e-11,2.63

cd /data/disk01/home/jasonfan/mapscript
./plotMercator.py $3$outname  -c 0 -L Significance -m 0 -M 30 -o $3$pre$png --origin 83.6332 22.0145 4 4|grep Max >> cut.txt