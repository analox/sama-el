#include <stdlib.h>
#include <stdio.h>
#include <stdarg.h>
#include <fcntl.h>

extern void clcdgrad_(float *, 
					   float *,
					   float *,
					   float *,
					   float *,
					   float *);

extern void shape_(float *, char *);

extern void eul2dplt_(float *, char *);

/*
How to Call the ShapePress Program

./ShapePress filename POPSIZE


*/


int main(int argc,char *argv[])
{
  int ND = 24;
  int POPSIZEMAX = 60;
  float alphIn[ND];
  float clOut;
  float cdOut;
  float gOut[ND];
  float rsdEu;
  float rsdAdj;
  char dname[8];
  FILE *fmax;
  FILE *fmin;
  FILE *ftyp;
  FILE *fpop;
  FILE *fduh;
  //char string[2000];
  float alphMax[ND];
  float alphMin[ND];
  float alph000[ND];
  float pop[POPSIZEMAX][ND+1];
  float printpop[ND+10];
  //char *token;
  char populationfit[100];
  char fout[100];
  int i,j,k;
  int POPSIZE; 

  sprintf(populationfit, "%s",argv[1]);
  printf("populationfit: %s\n",argv[1]);
  
  POPSIZE=atoi(argv[2]);
  
  for (i=0; i<ND; i++) alph000[i]=0;

  fmax = fopen("minmax24.dat","r");
  ftyp = fopen("ex-alpha24.dat","r");
  
  fpop = fopen(populationfit,"r");
  
  for (i=0; i<ND; i++) {
	fscanf(fmax," %f %f",&alphMin[i],&alphMax[i]);
	fscanf(ftyp," %f",&alphIn[i]);
  }
  
  
  /*
  while ((fgets (string, 2000, fpop)) != NULL)
  {
  	//printf ("<%s>\n", string);
  	j = 0;
    	i = 0;
    	token = strtok (string, " ");
    	//printf ("token <%s>\n", token);
 	do {
 		if (i > ND)
 		{
 			j++;
 			break;
 		}
 		
 		pop[j][i] = atof(token);
		i++;
    					
	} while ((token = strtok ((char *)NULL, " ")) != (char *)NULL);
  }
  
  for (j=0; j<POPSIZE; j++) {
  	printf(" %f %f %f %f %f %f %f %f %f %f
	%f %f %f %f %f %f %f %f %f %f
	%f %f %f %f %f %f %f %f %f %f
	%f %f %f %f %f %f %f %f %f %f
	%f %f %f %f %f %f %f %f %f %f
	%f %f %f %f %f %f %f %f %f %f
	%f\n",
	pop[j][0],
	pop[j][1], pop[j][2], pop[j][3], pop[j][4], pop[j][5], 
	pop[j][6], pop[j][7], pop[j][8], pop[j][9], pop[j][10],
	pop[j][11], pop[j][12], pop[j][13], pop[j][14], pop[j][15],
	pop[j][16], pop[j][17], pop[j][18], pop[j][19], pop[j][20],
	pop[j][21], pop[j][22], pop[j][23], pop[j][24], pop[j][25],
	pop[j][26], pop[j][27], pop[j][28], pop[j][29], pop[j][30],
	pop[j][31], pop[j][32], pop[j][33], pop[j][34], pop[j][35],
	pop[j][36], pop[j][37], pop[j][38], pop[j][39], pop[j][40],
	pop[j][41], pop[j][42], pop[j][43], pop[j][44], pop[j][45],
	pop[j][46], pop[j][47], pop[j][48], pop[j][49], pop[j][50],
	pop[j][51], pop[j][52], pop[j][53], pop[j][54], pop[j][55],
	pop[j][56], pop[j][57], pop[j][58], pop[j][59], pop[j][60]);
}
	
  */
  
  /*
  for (j=0; j<POPSIZE; j++) {
	fscanf(fpop,
	"%f %f %f %f %f %f %f %f %f %f
	%f %f %f %f %f %f %f %f %f %f
	%f %f %f %f %f %f %f %f %f %f
	%f %f %f %f %f %f %f %f %f %f
	%f %f %f %f %f %f %f %f %f %f
	%f %f %f %f %f %f %f %f %f %f
	%f",
	&pop[j][0],
	&pop[j][1], &pop[j][2], &pop[j][3], &pop[j][4], &pop[j][5], 
	&pop[j][6], &pop[j][7], &pop[j][8], &pop[j][9], &pop[j][10],
	&pop[j][11], &pop[j][12], &pop[j][13], &pop[j][14], &pop[j][15],
	&pop[j][16], &pop[j][17], &pop[j][18], &pop[j][19], &pop[j][20],
	&pop[j][21], &pop[j][22], &pop[j][23], &pop[j][24], &pop[j][25],
	&pop[j][26], &pop[j][27], &pop[j][28], &pop[j][29], &pop[j][30],
	&pop[j][31], &pop[j][32], &pop[j][33], &pop[j][34], &pop[j][35],
	&pop[j][36], &pop[j][37], &pop[j][38], &pop[j][39], &pop[j][40],
	&pop[j][41], &pop[j][42], &pop[j][43], &pop[j][44], &pop[j][45],
	&pop[j][46], &pop[j][47], &pop[j][48], &pop[j][49], &pop[j][50],
	&pop[j][51], &pop[j][52], &pop[j][53], &pop[j][54], &pop[j][55],
	&pop[j][56], &pop[j][57], &pop[j][58], &pop[j][59], &pop[j][60]);
	
	
	sprintf(fout, "shapePop%dChro%d",atoi(argv[1]),);
	shape_(alphIn,"shapenew");
	
  	printf(" %f %f %f %f %f %f %f %f %f %f
	%f %f %f %f %f %f %f %f %f %f
	%f %f %f %f %f %f %f %f %f %f
	%f %f %f %f %f %f %f %f %f %f
	%f %f %f %f %f %f %f %f %f %f
	%f %f %f %f %f %f %f %f %f %f
	%f\n",
	pop[j][0],
	pop[j][1], pop[j][2], pop[j][3], pop[j][4], pop[j][5], 
	pop[j][6], pop[j][7], pop[j][8], pop[j][9], pop[j][10],
	pop[j][11], pop[j][12], pop[j][13], pop[j][14], pop[j][15],
	pop[j][16], pop[j][17], pop[j][18], pop[j][19], pop[j][20],
	pop[j][21], pop[j][22], pop[j][23], pop[j][24], pop[j][25],
	pop[j][26], pop[j][27], pop[j][28], pop[j][29], pop[j][30],
	pop[j][31], pop[j][32], pop[j][33], pop[j][34], pop[j][35],
	pop[j][36], pop[j][37], pop[j][38], pop[j][39], pop[j][40],
	pop[j][41], pop[j][42], pop[j][43], pop[j][44], pop[j][45],
	pop[j][46], pop[j][47], pop[j][48], pop[j][49], pop[j][50],
	pop[j][51], pop[j][52], pop[j][53], pop[j][54], pop[j][55],
	pop[j][56], pop[j][57], pop[j][58], pop[j][59], pop[j][60]);
	
  }
  */
  
  for (j=0; j<POPSIZE; j++) {
	fscanf(fpop,
	"%f %f %f %f %f %f %f %f %f %f
	%f %f %f %f %f %f %f %f %f %f
	%f %f %f %f %f %f",
	&printpop[0],
	&printpop[1], &printpop[2], &printpop[3], &printpop[4], &printpop[5], 
	&printpop[6], &printpop[7], &printpop[8], &printpop[9], &printpop[10],
	&printpop[11], &printpop[12], &printpop[13], &printpop[14], &printpop[15],
	&printpop[16], &printpop[17], &printpop[18], &printpop[19], &printpop[20],
	&printpop[21], &printpop[22], &printpop[23], &printpop[24], &printpop[25] );
	
	/*
	for (k=0; k<ND; k++) {
		printpop[k] = printpop[k]/10e5;
	}
	*/
	
	printf("ppD%dEv%d\n",j,(int)printpop[24]);
	sprintf(fout, "ppD%dEvalu\n",j);
	printf("fout is %s\n", fout);
	shape_(printpop,fout);	
  }
  
  fcloseall();

  /* Example of how to call clcdgrad_ */
  /* Calculate for shape defined by alpha: */
  /*   float clOut = lift coeff */
  /*   float cdOut = drag coeff */
  /*   float gOut[ND]  = drag gradient */
  /*   float rsdEu = Euler residual */
  /*   float rsdAdj = Adjoint residual */
  
  // clcdgrad_(alphIn,&clOut,&cdOut,gOut,&rsdEu,&rsdAdj);

  /* Example of how to call shape_     */
  /* Generate shape coordinates defined by alpha. */
  /* Input:                            */
  /*    Argument 1 = float alpha[ND]   */
  /*    Argument 2 = character name[8] */
  /* Output files:                     */
  /*    <name>.plt                     */
  shape_(alph000,"shape000");
  shape_(alphMax,"shapemax");
  shape_(alphMin,"shapemin");
  shape_(alphIn,"shapenew");
  
  /* Example of how to call eul2dplt_ */
  /* Solve the Euler field solution for shape defined by alpha. */
  /* Input:                            */
  /*    Argument 1 = float alpha[ND]   */
  /*    Argument 2 = character name[8] */
  /* Output files:                     */
  /*  	<name>-cpl.plt                 */
  /*  	<name>-cpu.plt                 */
  /*  	<name>-fld.plt                 */

  //eul2dplt_(alphMax,"example");

  return 0; 
}


