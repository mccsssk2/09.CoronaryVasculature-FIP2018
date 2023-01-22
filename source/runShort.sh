#!/bin/bash

	make veryclean
	make vasc

for ((b=1; b <=100 ; b++))
do
#	rm -rf raw*
	./vasc 4 6
	mv rawSegsRCA.data rawSegsRCA.data.$b
	mv rawSegsLA.data  rawSegsLA.data.$b
	mv segmentsRadius3DRCA.vtk segmentsRadius3DRCA$b.vtk
	mv segmentsRadius3DLA.vtk segmentsRadius3DLA$b.vtk
done
