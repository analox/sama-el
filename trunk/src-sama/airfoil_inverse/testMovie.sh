#!/bin/sh

for i in `seq 1 15`;
do
	echo "running number $i"
	cp ./movie/shape$i/design$i.dat ex-alpha24.dat
	./lana 0.5 2.0
	cp result.dat ./movie/shape$i/
done
