#!/bin/bash

make veryclean
make veryclean
make vasc

rm -rf RCA
rm -rf LAD
rm -rf LCX

mkdir -p RCA
mkdir -p LAD
mkdir -p LCX

cp vasc RCA/
cp vasc LAD/
cp vasc LCX/

cp runShortRCA.sh RCA/
cp runShortLAD.sh LAD/
cp runShortLCX.sh LCX/

cd RCA
./runShortRCA.sh &> outt.datata &
cd ..
cd LAD
./runShortLAD.sh &> outt.datata &
cd ..
cd LCX
./runShortLCX.sh &> outt.datata &
cd ..

