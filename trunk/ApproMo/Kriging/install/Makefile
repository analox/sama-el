#  makefile last modified 6/1/00
#
#  Fill in the compiler you use.  We used g++.
CC=g++
#  Fill in the flags you wish for that compiler
CCFLAGS=-g -Wall -O
#These are for the CLAPACK, BLAS, and F2C libraries.  You should modify them
#to match whatever you have named these libraries.
LIBS=-L./ -lLAPACK -lBLAS -lf2c

# To compile the krigifier with some other main program than testkrigify.cc, simply go down to the\
# entry for testkrigify, and replace all occurences of testkrigify.cc with your own .cc filename\
# feel free to change this Makefile to suit your needs

default: demo1

demos: demo1 demo2 demo3 demo4 demo5 demo6 plotpoints

kriggui: kriggui.o rngs.o rvgs.o krigify.o gamma.o kriggui.tcl plotpoints
	$(CC) $(CCFLAGS) kriggui.o rngs.o rvgs.o krigify.o\
	gamma.o $(LIBS) -o kriggui

demo1: krigify.o demo1.o rngs.o rvgs.o krig.o\
	gamma.o ParamEstimate.o DirectSearch.o\
	PatternSearch.o CompassSearch.o approx.o objective.o
	$(CC) $(CCFLAGS) demo1.o krig.o krigify.o rngs.o rvgs.o\
	gamma.o ParamEstimate.o PatternSearch.o CompassSearch.o\
	approx.o DirectSearch.o objective.o $(LIBS) -o demo1

demo2: krigify.o demo2.o rngs.o rvgs.o krig.o\
	gamma.o ParamEstimate.o DirectSearch.o\
	PatternSearch.o CompassSearch.o approx.o objective.o
	$(CC) $(CCFLAGS) demo2.o krig.o krigify.o rngs.o rvgs.o\
	gamma.o ParamEstimate.o PatternSearch.o CompassSearch.o\
	approx.o DirectSearch.o objective.o $(LIBS) -o demo2

demo3: krigify.o demo3.o rngs.o rvgs.o gamma.o
	$(CC) $(CCFLAGS) demo3.o krigify.o rngs.o rvgs.o gamma.o\
	$(LIBS) -o demo3

demo4: krigify.o demo4.o rngs.o rvgs.o krig.o\
	gamma.o ParamEstimate.o DirectSearch.o\
	PatternSearch.o CompassSearch.o approx.o objective.o
	$(CC) $(CCFLAGS) demo4.o krig.o krigify.o rngs.o rvgs.o\
	gamma.o ParamEstimate.o PatternSearch.o CompassSearch.o\
	approx.o DirectSearch.o objective.o $(LIBS) -o demo4

demo5: krigify.o demo5.o rngs.o rvgs.o krig.o\
	gamma.o ParamEstimate.o gls.o DirectSearch.o\
	PatternSearch.o CompassSearch.o approx.o objective.o
	$(CC) $(CCFLAGS) demo5.o krig.o krigify.o rngs.o rvgs.o\
	gamma.o ParamEstimate.o PatternSearch.o CompassSearch.o\
	approx.o gls.o DirectSearch.o objective.o $(LIBS) -o demo5

demo6: gls.o demo6.o rngs.o rvgs.o gamma.o krig.o ParamEstimate.o\
	PatternSearch.o CompassSearch.o approx.o DirectSearch.o objective.o
	$(CC) $(CCFLAGS) demo6.o gls.o rngs.o rvgs.o krig.o gamma.o\
	ParamEstimate.o PatternSearch.o CompassSearch.o approx.o\
	DirectSearch.o objective.o $(LIBS) -o demo6

demoSVD: krigifySVD.o demo.o rngs.o rvgs.o krigSVD.o\
	gamma.o ParamEstimate.o DirectSearch.o\
	PatternSearch.o CompassSearch.o approx.o
	$(CC) $(CCFLAGS) demo.o krigSVD.o krigifySVD.o rngs.o rvgs.o\
	gamma.o ParamEstimate.o PatternSearch.o CompassSearch.o\
	approx.o DirectSearch.o $(LIBS) -o demoC

testkrigify: testkrigify.cc krigify.o gamma.o rngs.o rvgs.o
	$(CC) $(CCFLAGS) testkrigify.cc krigify.o gamma.o rngs.o\
	rvgs.o $(LIBS) -o testkrigify 

krigify.o: krigify.cc krigify.h maps_general.h vec.h cppmat.h rngs.h\
	rvgs.h  gamma.h Dyn_alloc.h
	$(CC) $(CCFLAGS) -c krigify.cc

#This for the pre-Cholesky version of the krigifier
krigifySVD.o: krigifySVD.cc krigify.h maps_general.h vec.h cppmat.h rngs.h\
	rvgs.h  gamma.h
	$(CC) $(CCFLAGS) -c krigifySVD.cc

rngs.o: rngs.c rngs.h
	$(CC) $(CCFLAGS) -c rngs.c

rvgs.o: rvgs.c rvgs.h
	$(CC) $(CCFLAGS) -c rvgs.c

rbf.o: rbf.cc rbf.h gamma.h approx.h
	$(CC) $(CCFLAGS) -c rbf.cc

gamma.o: gamma.cc gamma.h
	$(CC) $(CCFLAGS) -c gamma.cc

krig.o: krig.cc krig.h ParamEstimate.h gamma.h CompassSearch.h approx.h Dyn_alloc.h
	$(CC) $(CCFLAGS) -c krig.cc

gls.o: gls.cc gls.h gamma.h krig.h approx.h ParamEstimate.h CompassSearch.h Dyn_alloc.h
	$(CC) $(CCFLAGS) -c gls.cc

#This is for the pre-Cholesky version of the kriging approximation
krigSVD.o: krigSVD.cc krig.h ParamEstimate.h gamma.h CompassSearch.h approx.h
	$(CC) $(CCFLAGS) -c krigSVD.cc

ParamEstimate.o: ParamEstimate.h ParamEstimate.cc PatternSearch.h maps_general.h
	$(CC) $(CCFLAGS) -c ParamEstimate.cc

PatternSearch.o: maps_general.h vec.h cppmat.h  PatternSearch.h PatternSearch.cc
	$(CC) $(CCFLAGS) -c PatternSearch.cc

CompassSearch.o: maps_general.h vec.h cppmat.h  CompassSearch.h CompassSearch.cc
	$(CC) $(CCFLAGS) -c CompassSearch.cc

DirectSearch.o: maps_general.h vec.h cppmat.h  DirectSearch.h\
	DirectSearch.cc objective.h Dyn_alloc.h
	$(CC) $(CCFLAGS) -c DirectSearch.cc

approx.o: approx.cc maps_general.h approx.h
	$(CC) $(CCFLAGS) -c approx.cc

demo1.o: demo1.cc krig.h krigify.h rngs.h rvgs.h
	$(CC) $(CCFLAGS) -c demo1.cc

demo2.o: demo2.cc krig.h krigify.h rngs.h rvgs.h
	$(CC) $(CCFLAGS) -c demo2.cc

demo3.o: demo3.cc krigify.h rngs.h rvgs.h gamma.h
	$(CC) $(CCFLAGS) -c demo3.cc

demo4.o: demo4.cc krig.h krigify.h rngs.h rvgs.h
	$(CC) $(CCFLAGS) -c demo4.cc

demo5.o: demo5.cc gls.h krig.h krigify.h rngs.h rvgs.h
	$(CC) $(CCFLAGS) -c demo5.cc

demo6.o: demo6.cc krig.h gls.h rngs.h rvgs.h
	$(CC) $(CCFLAGS) -c demo6.cc

plotpoints: plotpoints.c
	gcc -O plotpoints.c -o plotpoints

kriggui.o: krigify.h rngs.h rvgs.h kriggui.cc
	$(CC) $(CCFLAGS) -c kriggui.cc

clean: 
	rm -f *.o || /bin/true
	rm -f *~ || /bin/true
	rm -f core || /bin/true

