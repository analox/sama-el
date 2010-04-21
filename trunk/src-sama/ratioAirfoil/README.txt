To compile: make ratioairfoilObjIn

Example of usage: see file ratioairfoilObjIn.c

Input: 24 design variables from file ex-alpha24.dat, the lower and upper bound
of these 24 design variables are as written in minmax24.dat. In addition to
these variables, two environmental parameters, i.e. mach number and angle of
attack also can be input.

To run : ./ratioairfoilObjIn <mach_number> <angle_of_attack>
e.g.: ./ratioairfoilObjIn 0.5 2.0
