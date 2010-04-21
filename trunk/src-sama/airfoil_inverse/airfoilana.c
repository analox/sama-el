#include <stdlib.h>
#include <stdio.h>
#include <stdarg.h>


extern void pdgradin_(float *,
                                           float *,
                                           float *,
                                           float *,
                                           float *,
                                           float *,
                                           float *);

// mach = static variable with noise
// aoa = 0.0, no noise
// alphaIn = design variables, no noise

void airfoilana(float *mach, float *aoa, float *alphIn, float *pdOut, float *gOut, float * rsdEu, float *rsdAdj, int ND)
{
  int i;
  pdgradin_(mach,aoa,alphIn,pdOut,gOut,rsdEu,rsdAdj);

}

