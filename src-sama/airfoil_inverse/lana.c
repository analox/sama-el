#include <stdlib.h>
#include <stdio.h>
#include <stdarg.h>
#include <time.h>


int main(int argc, char**argv)
{
  time_t date1,date2;
  int ND = 24;
  int i;

  float alphIn[ND];
  float pdOut;
  float gOut[ND];
  float rsdEu;
  float rsdAdj;

  float mach,aoa;

  float alphMax[ND];
  float alphMin[ND];
  float alph000[ND];
  float alph001[ND];
  float test[ND];

  FILE *ftyp;
  //FILE *fmax;
  FILE *resultOut;

  ftyp = fopen("ex-alpha24.dat","r");
  //fmax = fopen("minmax24.dat","r");
  fprintf(stderr,"FILE is OPEN!!!\n");

  //aoa = 2.0;
  //mach =0.5;

  if(argc>1) {
  	mach = atof(argv[1]);
  	if(argc>2) aoa = atof(argv[2]);
  }

  // load some alpha values
  for (i=0; i<ND; i++) {
	//fscanf(fmax," %f %f",&alphMin[i],&alphMax[i]);
	fscanf(ftyp," %f",&alphIn[i]);
	//alphIn[i]/=1000000;
	//alph000[i]=0;
	//alph001[i]=1/1024;
  }


  printf("Evaluating alphIn for mach=%f, aoa=%f........\n", mach, aoa);
  airfoilana(&mach, &aoa, alphIn, &pdOut, gOut, &rsdEu, &rsdAdj, ND);
  printf("\npdOut=%f\n", pdOut);

  resultOut = fopen("result.dat","w");
  fprintf(resultOut,"%f\n", pdOut);
  fclose(resultOut);

  return 0;
}


