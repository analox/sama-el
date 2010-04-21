#!/bin/sh

./run.sh $1 params/param.rastrigin.dat > rastrigin.txt
mv logs Rastrigin
mkdir logs

./run.sh $1 params/param.ackley.dat > ackley.txt
mv logs Ackley
mkdir logs

./run.sh $1 params/param.griewank.dat > griewank.txt
mv logs Griewank
mkdir logs

./run.sh $1 params/param.rosenbrock.dat > rosenbrock.txt
mv logs Rosenbrock
mkdir logs

./run.sh $1 params/param.weierstrass.dat > weierstrass.txt
mv logs Weierstrass
mkdir logs

./run.sh $1 params/param.grierose.dat > grierose.txt
mv logs Grie+Rosen
mkdir logs

