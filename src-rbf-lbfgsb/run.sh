#!/bin/sh

rm logs/*
for i in `seq 1 $1`;
do
	./sama
	cp logs/bestFile.dat logs/bestFile-$i.dat
	cp logs/bestVecFile.dat logs/bestVecFile-$i.dat
done
rm logs/bestFile.dat logs/bestVecFile.dat

