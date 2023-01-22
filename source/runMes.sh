#!/bin/bash


make mes

mkdir -p distances
cp mes distances
cd distances
./mes 0.5 &> 0.5.data &
./mes 1.0 &> 1.0.data &
./mes 2.0 &> 2.0.data &
./mes 3.0 &> 3.0.data &
./mes 4.0 &> 4.0.data &
