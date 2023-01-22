#!/bin/bash

LIMIT=10

touch viz
make viz

for ((a=1; a <= LIMIT ; a++))
do
	cp viz ORCA_Oct4/ALL_$a/
	cd ORCA_Oct4/ALL_$a
	./viz
	cd ../..
done

