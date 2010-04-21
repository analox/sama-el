#!/bin/sh

rm logs/*
for i in `seq 1 $1`;
do
	./sama-airfoil $2
	cp logs/bestFile.dat logs/bestFile-$i.dat
	cp logs/freqFile.dat logs/freqFile-$i.dat
	cp logs/bestVecFile.dat logs/bestVecFile-$i.dat
done
rm logs/bestFile.dat logs/bestVecFile.dat logs/freqFile.dat

