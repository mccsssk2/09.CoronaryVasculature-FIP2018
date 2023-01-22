#!/bin/bash

for ((b=1; b <=100 ; b++))
do

./vasc 1 

mv meandParent.txt meandParent$b.txt

done
