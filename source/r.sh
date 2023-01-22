#!/bin/bash
LIMIT=10

for ((a=1; a <= LIMIT ; a++))
do
	DIR=RCA
	cp ${DIR}_$a/segmentsRadius3D.vtk segmentsRadius3D_${DIR}_$a.vtk
done


for ((a=1; a <= LIMIT ; a++))
do
	DIR=LAD
	cp ${DIR}_$a/segmentsRadius3D.vtk segmentsRadius3D_${DIR}_$a.vtk
done

for ((a=1; a <= LIMIT ; a++))
do
	DIR=LCX
	cp ${DIR}_$a/segmentsRadius3D.vtk segmentsRadius3D_${DIR}_$a.vtk
done

for ((a=1; a <= LIMIT ; a++))
do
	DIR=ALL
	cp ${DIR}_$a/segmentsRadius3D.vtk segmentsRadius3D_${DIR}_$a.vtk
done

