#!/bin/bash

for i in `seq 1 50`;
do
	cp viz ALL_$i
	cd ALL_$i
	./viz
	cd ..
done
