#!/bin/bash

cmd1="mkdir validation";
cmd2="mkdir validation/manual";
cmd3="mkdir validation/raw";

echo $cmd1
eval $cmd1
echo $cmd2
eval $cmd2
echo $cmd3
eval $cmd3

for i in 1 2 3 4 5 6 7 8 9 10
do

var1=$(echo "scale=2; $i/100" | bc)
var2=$(echo "scale=2; 1.1-$i/10" | bc)

cmd4="./cvtGenerator $((100*$i)) validation/manual/$i.mha validation/raw/$i.mha $var1 $var2";

echo $cmd4
eval $cmd4

done


