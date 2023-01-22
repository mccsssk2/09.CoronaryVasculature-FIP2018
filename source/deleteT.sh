#!/bin/bash


ulimit -s unlimited

start=892909
stop=892958

for (( j=${start}; j<=${stop}; j++ ));
do
	scancel $j
done
