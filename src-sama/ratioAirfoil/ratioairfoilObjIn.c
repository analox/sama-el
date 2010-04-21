#include <stdlib.h>
#include <stdio.h>
#include <stdarg.h>
#include <time.h>

extern void ratiogradin_(float *, //mach
					   float *, //aoa
					   float *, //alphIn
					   float *, //clOut
					   float *, //cdOut
					   float *, //gRat
					   float *, //rsdEu
					   float *, //rsdAdjCl
					   float *);//rsdAdjCd

//extern void shape_(float *, char *);

//extern void eul2dplt_(float *, char *, float *, float *);

void ratioairfoilObjIn();

int main(int argc , char** argv)
{
  int ND = 24;
  float alphIn[ND];
  float clOut, cdOut;
  float gRat[ND];
  float rsdEu;
  float rsdAdjCl, rsdAdjCd;
  char dname[8];
  FILE *fmax;
  FILE *fmin;
  FILE *ftyp;
  FILE *fout;
  float alphMax[ND];
  float alphMin[ND];
  float alph000[ND];
  float mach,aoa;
  int i,j;
  time_t date1 , date2;

  srand((unsigned int) time(NULL));

  for (i=0; i<ND; i++) alph000[i]=0;

  //fmax = fopen("minmax24.dat","r");
  ftyp = fopen("ex-alpha24.dat","r");

  /* load some alpha values */
  for (i=0; i<ND; i++) {
	//fscanf(fmax," %f %f",&alphMin[i],&alphMax[i]);
	fscanf(ftyp," %f",&alphIn[i]);
  }

  fcloseall();

  /* Example of how to call ratiograd_ */
  /* Calculate for shape defined by alpha: */
  /*   float clOut = lift coefficeint */
  /*   float cdOut = drag coefficeint */
  /*   float gRat[ND] = drag-over-lift ratio gradient */
  /*   float rsdEu = Euler residual */
  /*   float rsdAdjCl = Adjoint residual for lift */
  /*   float rsdAdjCl = Adjoint residual for drag */

  mach=0.5;
  aoa=2.0;
  if(argc>1) {
	mach=atof(argv[1]);
	aoa=atof(argv[2]);
  }

  printf("start calculating for mach = %f , aoa = %f \n",mach,aoa);
  time(&date1);
  ratioairfoilObjIn(&mach,&aoa,alphIn,&clOut,&cdOut,gRat,&rsdEu,&rsdAdjCl,&rsdAdjCd,ND);
  //printf("\tclOut\t= %f\tcdOut\t= %f\n",clOut,cdOut);
  //printf("\trsdEu\t= %f\n",rsdEu);
  //printf("\trsdAdjCl\t= %f\trsdAdjCl\t= %f\n\n",rsdAdjCl,rsdAdjCd);

  /* Example of how to call shape_     */
  /* Generate shape coordinates defined by alpha. */
  /* Input:                            */
  /*    Argument 1 = float alpha[ND]   */
  /*    Argument 2 = character name[8] */
  /* Output files:                     */
  /*    <name>.plt                     */
  //shape_(alphIn,"shapenew");

  /* Example of how to call eul2dplt_ */
  /* Solve the Euler field solution for shape defined by alpha. */
  /* Input:                            */
  /*    Argument 1 = float alpha[ND]   */
  /*    Argument 2 = character name[8] */
  /* Output files:                     */
  /*  	<name>-cpl.plt                 */
  /*  	<name>-cpu.plt                 */
  /*  	<name>-fld.plt                 */

  //eul2dplt_(alphIn,"example",&mach,&aoa);

  time(&date2);
  printf("Elapsed time: %d\n",(int)date2-(int)date1);

  fout=fopen("pdResult.dat","w+");
  fprintf(fout,"%f %f %f",cdOut, clOut,cdOut/clOut);

  return 0;
}


void ratioairfoilObjIn(float *mach, float *aoa, float *alphIn, float *clOut, float *cdOut, float *gRat, float *rsdEu, float *rsdAdjCl, float *rsdAdjCd, int ND)
{
  int i;

  ratiogradin_(mach,aoa,alphIn,clOut,cdOut,gRat,rsdEu,rsdAdjCl,rsdAdjCd);

  printf("mach=%f, aoa=%f, CL=%f, CD=%f, RSDEU=%f, RSDADJCL=%f, RSDADJCD=%f\n",mach[0],aoa[0],clOut[0],cdOut[0],rsdEu[0],rsdAdjCl[0],rsdAdjCd[0]);
  printf("Ratio Gradient = \n");
  for (i=0; i<ND; i++) {
  	printf("         %f\n",gRat[i]);
  }
}


