#!/bin/bash

for ((b=1; b <=100 ; b++))
do

./vasc 3

mv meandParent.txt meandParent$b.txt

done
