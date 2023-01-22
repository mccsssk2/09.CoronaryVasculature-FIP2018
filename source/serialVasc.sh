#!/bin/bash

# As suggested by Sharcnet, you make your directories, and submit individual jobs.

ulimit -s unlimited

# prep the data.
# you cannot do division easily in bash. so this is an integer that goes from 0 to 10 representing 0% to 100% coverage of area by ischemia patches.
# my job array is actually several separate jobs.
# module unload intel/12.1.3
# module unload openmpi/intel/1.6.2
# module load gcc/6.3.0


declare -a array=("RCA" "LAD" "LCX" "ALL")


# get length of an array
arraylength=${#array[@]}

for (( j=4; j<${arraylength}+1; j++ ));
do

# use for loop to read all values and indexes
for (( i = 1 ; i <=50; i++ ));
	  do 
		mkdir -p ${array[$j-1]}_$i
#		sqsub -r 48h -q serial --mpp=16G -o ${array[$j-1]}_$i/out.$i$j ./vasc 4 6 
		sbatch -o ${array[$j-1]}_$i/out.$i$j submit.sh
	done
done
