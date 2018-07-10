#!/bin/sh
#export PATH="~/miniconda2/bin:$PATH";
export LD_LIBRARY_PATH="~/miniconda2/lib:$LD_LIBRARY_PATH"
export PATH="~/miniconda2/envs/myconda/bin:$PATH"
eval `$HOME/utils/ape/current/ape sh externals`
eval `~/aerie/install/bin/hawc-config --env-sh`
export LD_LIBRARY_PATH=${CONDA_PREFIX}/lib:/data/disk01/software/hawc/ApeInstalled/x86_64/External/gsl/1.16/lib/:${LD_LIBRARY_PATH}
export CONFIG_HAWC=/data/disk01/home/jasonfan/config-hawc
#cd $1
for j in {12..13}
do
for i in `find $1$j -name "Aft*"|cut -f 5 -d _|sort|uniq`
do
echo $i
xcdf select "1" -o $1$j"/"$i"_"$j".xcd" `find $1$j -name "Aft*$i*"|sort`
echo "finish "$i
rm -rf `find $1$j -name "Aft*$i*"`
done
done

