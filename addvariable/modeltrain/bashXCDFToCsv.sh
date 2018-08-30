#!/bin/bash
#SBATCH --mem=200
#SBATCH -t 0-20:00:00

eval `$HAWCSOFT/setup.sh`
export PATH="~/miniconda2/bin:$PATH" 
export LD_LIBRARY_PATH="~/miniconda2/lib:$LD_LIBRARY_PATH"
source activate myconda
eval `$HOME/utils/ape/current/ape sh externals` 
eval `~/aerie/install/bin/hawc-config --env-sh`
export CONFIG_HAWC=/data/disk01/home/tcapistran/software/config-hawc
source activate
export LD_LIBRARY_PATH=${CONDA_PREFIX}/lib:/data/disk01/software/hawc/ApeInstalled/x86_64/External/gsl/1.16/lib/:${LD_LIBRARY_PATH}



ALLVARIABLES="0"
AllHADRON="1"

VERSION=""

DIROUT="/data/disk01/home/jasonfan"$VERSION
mkdir -p $DIROUT

PATHGAMMAFILE="/data/scratch/userspace/pretz/daqsim-reconstruction/output/daqsim-baseline-take4/"
GAMMAFILE=$PATHGAMMAFILE"gamma.xcd"
#CRFILESMC="carbon.xcd helium.xcd iron.xcd magnesium.xcd neon.xcd oxygen.xcd proton.xcd silicon.xcd"
CRFILESMC=""
PATHHADRONFILE="/lustre/scr01/userspace/pretz/run005214-sample-full/"
HADRONFILE=$PATHHADRONFILE"run005214.dec20.xcd"
#HADRONFILE=""

INFILE=$GAMMAFILE
if [ $AllHADRON == "1" ]
then
    for line in $CRFILESMC
    do
        INFILE=$INFILE" "$PATHGAMMAFILE"/"$line
    done
fi
CUT="(rec.angleFitStatus==0)&&(rec.coreFitStatus==0)&&(rec.nChAvail >= 700)&&(rec.coreFiduScale <= 100)&&(rec.zenithAngle < 0.785)"

VARIABLESRD="rec.logGPEnergy,rec.nChAvail,rec.nHit,rec.nHitSP10,rec.nHitSP20,rec.nTankAvail,rec.nTankHit,rec.windowHits,rec.planeNDOF,rec.SFCFNDOF,rec.CxPE40XnCh,rec.coreFiduScale,sweets.SWgt,sweets.TWgt,rec.zenithAngle,rec.azimuthAngle,rec.planeChi2,rec.coreX,rec.coreY,rec.logCoreAmplitude,rec.coreFitUnc,rec.SFCFChi2,rec.logNNEnergyV2,rec.fAnnulusCharge0,rec.fAnnulusCharge1,rec.fAnnulusCharge2,rec.fAnnulusCharge3,rec.fAnnulusCharge4,rec.fAnnulusCharge5,rec.fAnnulusCharge6,rec.fAnnulusCharge7,rec.fAnnulusCharge8,rec.protonlheEnergy,rec.protonlheLLH,rec.gammalheEnergy,rec.gammalheLLH,rec.logMaxPE,rec.logNPE,rec.CxPE40,rec.CxPE40SPTime,rec.LDFAge,rec.LDFAmp,rec.LDFChi2,rec.PINC,rec.disMax"


TAG=`basename $HADRONFILE .xcd`
OUTPUTFILE=$DIROUT"/"$TAG".csv"
echo "INPUT: "$HADRONFILE
echo "VARIABLES: "$VARIABLESRD
if [ $ALLVARIABLES == "1" ]
then
    OUTPUTFILE=$DIROUT"/All"$TAG".csv"
    echo "OUTPUT: "$OUTPUTFILE
    xcdf select "$CUT" $HADRONFILE | xcdf csv > $OUTPUTFILE
else
    echo "OUTPUT: "$OUTPUTFILE
    xcdf select "$CUT" $HADRONFILE | xcdf select-fields "$VARIABLESRD" | xcdf csv > $OUTPUTFILE
fi


VARIABLESMC="sweets.oneWgt,sweets.IWgt,sweets.TWgt,sweets.BWgt,rec.nChAvail,rec.nHit,rec.nHitSP10,rec.nHitSP20,rec.nTankAvail,rec.nTankHit,rec.windowHits,rec.planeNDOF,rec.SFCFNDOF,rec.CxPE40XnCh,rec.coreFiduScale,mc.corsikaParticleId,rec.zenithAngle,rec.azimuthAngle,rec.dec,rec.ra,rec.planeChi2,rec.coreX,rec.coreY,rec.logCoreAmplitude,rec.coreFitUnc,rec.SFCFChi2,rec.logNNEnergyV2,rec.fAnnulusCharge0,rec.fAnnulusCharge1,rec.fAnnulusCharge2,rec.fAnnulusCharge3,rec.fAnnulusCharge4,rec.fAnnulusCharge5,rec.fAnnulusCharge6,rec.fAnnulusCharge7,rec.fAnnulusCharge8,rec.protonlheEnergy,rec.protonlheLLH,rec.gammalheEnergy,rec.gammalheLLH,rec.logMaxPE,rec.logNPE,rec.CxPE40,rec.CxPE40SPTime,rec.LDFAge,rec.LDFAmp,rec.LDFChi2,rec.PINC,rec.disMax,mc.delAngle,mc.logEnergy"
for line in $INFILE
do
    TAG=`basename $line .xcd`
    OUTPUTFILE=$DIROUT"/"$TAG".csv"
#    TAG=${line%.xcd}.csv
    echo "INPUT: "$line
    if [ $ALLVARIABLES == "1" ]
    then
        OUTPUTFILE=$DIROUT"/All"$TAG".csv"
        echo "OUTPUT: "$OUTPUTFILE
        xcdf select "$CUT" $line | xcdf csv > $OUTPUTFILE
    else
        echo "OUTPUT: "$OUTPUTFILE
        xcdf select "$CUT" $line | xcdf select-fields "$VARIABLESMC" | xcdf csv > $OUTPUTFILE
    fi
done


