#!/bin/bash

LIMIT=16

for ((a=11; a <= LIMIT ; a++))
do
	DIR=ALL
	mkdir -p ${DIR}_$a
	cp vasc ${DIR}_$a
	cd ${DIR}_$a
	./vasc 4 6 &> out.datat &
	cd ..
done



