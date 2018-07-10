temp=0
cat submap.sh | while read p
do
tee=" &"
ttt="nohup "
i=$ttt$p$tee
eval $i 
#echo $p &
temp=$((temp+1))
echo $temp
if (( $temp % 48 == 0 )); then wait; fi
done 
wait
