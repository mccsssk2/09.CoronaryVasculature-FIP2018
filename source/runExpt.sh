#!/bin/bash

for ((a=0; a <= 5 ; a++))
do
	./expt $a &> out.datat.$a &
	cd ..
done



