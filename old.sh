#!/bin/sh
: '
num_procs=$1
num_jobs="\j" 
cat temp.txt | while read p
do
while (( ${num_jobs@P} >= num_procs  )); do
    wait -n
done
  
tee=" &"
ttt="nohup "
i=$ttt$p$tee
eval $i 

done 
'

open_sem(){
    mkfifo pipe-$$
    exec 3<>pipe-$$
    rm pipe-$$
    local i=$1
    for((;i>0;i--)); do
        printf %s 000 >&3
    done
}
run_with_lock(){
    local x
    read -u 3 -n 3 x && ((0==x)) || exit $x
    (
    "$@" 
    printf '%.3d' $? >&3
    )&
}

N=8
open_sem $N
cat temp.txt | while read p
do
    run_with_lock $p
done 