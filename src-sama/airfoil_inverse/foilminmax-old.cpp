#include <stdio.h>
#include <iostream.h>
#include <iomanip.h>
#include <fstream.h>
#include <ga/ga.h>
#include <time.h>
#include <stdlib.h>
#include <math.h>
#include <cstring>
#include <signal.h>
#include "sort.h"

#define dimension 24
#define phi 3.141592654
#define MAX_Indiv 2000
#define MAX_ITER 3

float foilObjective(GAGenome &);
float findWorstPdOutExact(float* genome, float lRange, float uRange, float fx) ;
float findWorstPdOutApprox(float* genome, float lRange, float uRange, float fx) ;
void  RBFLocalSearch(int traindatasize, int no_of_var, int kerneltype, float kernelhyp, float Regparam, int nf, int ineqnf, int eqnf, int maxiter, float *sminbounds, float *smaxbounds, float (*traindata)[dimension+2], float *startdesign, int maximize);
void  sort();
int checkRepeatDesign(float machApprox, float* genomeApprox);

float inPoints[MAX_Indiv][dimension+1+2]={0};
// inPoints[][0] contains Mach
// inPoints[][1] ... inPoints[][24] contains Design Variables 1-24
// inPoints[][25] contains fitness, pdout
// inPoints[][26] contains ecluidean difference between the design vector to the design point of interest

int nOfGenExact=1;   //number of generation to use exact fitness function
float aoa=2.0;   //currently fixed to this value


int nOfBit=16;
int population=50,generation=1000;         //if use term=2, means termination using convergence, number of generation has no effect on termination
float bestfitness=100.0;

short term=2;
float alphaMax[dimension];
float alphaMin[dimension];
int currGeneration = 0;
int currCount = 0;    //count the number of points in database


int countChromo=0;

int main(int argc, char **argv)
{
	//mutation and crossover probability, will only be used in main, so just use local variable
	float pmut   = 0.1;
	float pcross = 0.9;
	int i,j;
	cout<<"Starting "<<endl;
  	GABin2DecPhenotype map;
  	srand((unsigned int) time(NULL));

	ifstream fminmax("minmax24.dat",ios::in);
	if(fminmax==NULL) {
		cerr<<"Error in opening minmax24.dat, hence programm will exit here!!!"<<endl;
		exit(1);
	}

	float tempVal[dimension]={0.0};

  	/* load range of alpha values */
  	for (i=0; i<dimension; i++) {
		fminmax>>alphaMin[i]>>alphaMax[i];
		map.add(nOfBit, alphaMin[i], alphaMax[i]);
		//sigma[i] = 0.01;
		//sigma[i] = 0.01 * (alphaMax[i] - alphaMin[i]);
		//sigma[i] = 0.05 * (alphaMax[i] - alphaMin[i]);
 	 }
	fminmax.close();


// Create the template genome using the phenotype map we just made.

  	GABin2DecGenome genome(map, foilObjective);


// Now create the GA using the genome and run it.  We'll use sigma truncation
// scaling so that we can handle negative objective scores.

  	GASimpleGA ga(genome);
  	GASigmaTruncationScaling scaling;
  	ga.populationSize(population);
  	ga.nGenerations(generation);
  	ga.pMutation(pmut);
  	ga.pCrossover(pcross);
  	ga.scaling(scaling);
	ga.maximize();

	if(term==2) {
		ga.pConvergence(0.9);
		ga.nConvergence(100);
		ga.terminator(GAGeneticAlgorithm::TerminateUponConvergence);
	}
  	ga.initialize(time(0));


	ofstream bestGenomeFile, generationFile, logFile  ;
	logFile.open("logFile.dat",ios::out|ios::trunc);
	logFile.close();



	time_t currentTime;
	bool newBest;

	ofstream outFile,outFile1,outFile2,outFile3,outFile4;     //contains all chromosomes in every generation
	outFile.open("chromofile.dat",ios::trunc);
	outFile.close();
	
	
	while(!ga.done())
	{
		currGeneration=ga.generation();

		//chromofile.dat logs all chromosomes in every generation
		outFile.open("chromofile.dat",ios::out|ios::app);
		outFile<<ga.generation()<<"\t";
		float temp[dimension]={0};
		for(i=0; i<ga.population().size(); i++){
			outFile<<"chromo no. "<<i<<": ";
			genome=ga.population().individual(i);
			for(j=0;j<dimension;j++) {
				outFile<<setprecision(6)<<genome.phenotype(j)<< "\t";
				temp[j]+=genome.phenotype(j);
			}
			outFile<<endl;

		}
		outFile.close();

		//genFile.dat logs current generation number, for observation purpose only
		generationFile.open("genFile.dat",ios::out|ios::trunc);
		generationFile<<ga.generation();
		cout<<"Currently in generation "<<currGeneration<<endl;
		generationFile.close();

		//go to next generation
		ga.step();

		//get the best individual so far
		genome=ga.statistics().bestIndividual();

		//find whether it is the newer best
		newBest=false;
		for(i=0;i<dimension;i++) {
			if(tempVal[i]!=genome.phenotype(i)) {newBest=true; break;}

		}
		if(newBest) {
			time(&currentTime);
			logFile.open("logFile.dat",ios::out|ios::app);
			logFile<<asctime(localtime(&currentTime))<<endl;
			logFile<<"New Best individual found at generation "<<ga.generation()<<" :"<<endl;
			cout<<"New Best individual found at generation "<<ga.generation()<<" :"<<endl;
			bestGenomeFile.open("best-alpha24.dat",ios::out|ios::trunc);
			for(i=0;i<dimension;i++) {
				tempVal[i]=genome.phenotype(i);
				logFile<<genome.phenotype(i)<<endl;
				bestGenomeFile<<genome.phenotype(i)<<endl;
			}
			logFile<<endl;
			bestGenomeFile.close();
			logFile.close();
		}

		
	}

	cout<<"Finish!"<<endl;
	return 0;
}



/*short function for debugging purpose :)
float foilObjective(GAGenome & c) {
	cout<<"Evaluating"<<endl;
	return GARandomFloat();
}*/


float foilObjective(GAGenome& c) {
	
	countChromo=(countChromo++)%population;
	cout<<"chromo number : "<<countChromo<<"at generation :"<<currGeneration<<endl;
	
	
	int i,j;
	GABin2DecGenome & genome = (GABin2DecGenome &)c;
	float arr[dimension]={0.0};
	for(i=0;i<dimension;i++) {
		arr[i] = genome.phenotype(i);	
	}

	int fitVal = 0;
	float worstPdOut = 0;
	float errorPercentage = 300;    //assume the worst error accepted is 300%
	float midPoint = 0.5,lRange, uRange;

	float pdOut;
	//get f(x)
	char argu[40]="./lana 0.5 2.0 ";
	ofstream ftyp("ex-alpha24.dat",ios::out|ios::trunc);
	for( i=0;i<dimension;i++) {
		ftyp<<arr[i]<<'\n';
	}
	ftyp.close();
	system(argu);
	ifstream getPd("result.dat",ios::in);
	getPd>>pdOut;
	float fx = pdOut;
	getPd.close();

	int currCountHash;

	int repeatPoint = checkRepeatDesign(0.5, arr);
	if(!repeatPoint) {
		currCountHash=currCount%MAX_Indiv;
		inPoints[currCountHash][0] = 0.5;
		for(j=0; j<dimension; j++) {
			inPoints[currCountHash][j+1] = arr[j];
		}
		inPoints[currCountHash][dimension+1] = pdOut;
		currCount++;
	}

	float error=0.0;
	while( error < errorPercentage ) {
		fitVal+=10;
		lRange=midPoint-0.025;
		uRange=midPoint+0.025;
		if(currGeneration<nOfGenExact) worstPdOut = findWorstPdOutExact(arr, lRange, uRange,fx);
		else worstPdOut = findWorstPdOutApprox(arr, lRange, uRange,fx);
	
		error = (float)100.0*(worstPdOut-fx)/fx;
		if(error<0.0) error*=-1.0;
		
		
	}

	return fitVal;
}



float findWorstPdOutExact(float* genome, float lRange, float uRange, float fx) {
	cout<<"using exact"<<endl;
	int i, j;
	float pdOut;
	//get f(x)
//	char argu[40]="./lana 0.5 2.0 ";
	ofstream ftyp("ex-alpha24.dat",ios::out|ios::trunc);
//	for( i=0;i<dimension;i++) {
//		ftyp<<genome[i]<<'\n';
//	}
	ftyp.close();
//	system(argu);
	ifstream getPd("result.dat",ios::in);
//	getPd>>pdOut;
//	float fx = pdOut;
	getPd.close();

	

	float worstPd = 0;
	for(i=0; i<MAX_ITER;i++) {
		float mach2 = GARandomFloat(lRange, uRange);
		ftyp.open("ex-alpha24.dat",ios::out|ios::trunc);
		for( j=0;j<dimension;j++) {
			ftyp<<genome[j]<<'\n';
		}
		ftyp.close();
		char arg[40];
		sprintf(arg, "./lana %f 2.0", mach2);
		system(arg);
		getPd.open("result.dat",ios::in);
		getPd>>pdOut;
		getPd.close();

		if(pdOut>worstPd) worstPd = pdOut;

		//store in database
		int repeatPoint = checkRepeatDesign(mach2, genome);
  		if(!repeatPoint) {
  			int currCountHash=currCount%MAX_Indiv;
  			inPoints[currCountHash][0] = mach2;
  			for(j=0; j<dimension; j++) {
  				inPoints[currCountHash][j+1] = genome[j];
  			}
  			inPoints[currCountHash][dimension+1] = pdOut;
  			currCount++;
  		}
	}

	return worstPd;

}


float findWorstPdOutApprox(float* genome, float lRange, float uRange, float fx) {

	cout<<"approximating"<<endl;

	//use approximation here
	int i, j, k;
	float pdOut;

	/* Do not touch this, these are the parameters required by RBF Local Search */
	int kerneltype 	=1;
	float kernelhyp =1.0;
	float Regparam  =10E-10;
	int nf	   	=1;
	int ineqnf	=0;
	int eqnf 	=0;
	int maxiter	=200;
	// the worst case condition in a minimization problem for robust design is to do a local maximization
	int maximize    =1;
	int neighCount  =  100;        //maximum 500

	float infinityMin, infinityMax;
	float trust= 0.0;
	infinityMax = uRange;
	infinityMin = lRange;
	trust= (uRange-lRange) * 0.5;

	int currCountHash;
	float dataPoints[MAX_Indiv][dimension+2];

	float sMin[dimension+1], sMax[dimension+1];

	//get f(x)
//	char argu[40]="./lana 0.5 2.0 ";
	ofstream ftyp("ex-alpha24.dat",ios::out|ios::trunc);
//	for( j=0;j<dimension;j++) {
//		ftyp<<genome[j]<<'\n';
//	}
	ftyp.close();
//	system(argu);
	ifstream getPd("result.dat",ios::in);
//	getPd>>pdOut;
//	float fx = pdOut;
//	float diffsq;
	getPd.close();


	float faa, fxa, trustChange;
        int index;
        
        float genomeApprox[dimension];
	float machApprox;
	float resultApprox;     //faa
	float machInitial;
	float rho;
	float diffsq;
	float fax;
	
	for(k=0; k<MAX_ITER; k++) {
		cout<<"iter "<<k<<endl;
		int evalEuclidean = 0;
		for(i=0; i<dimension+1; i++) {
			if(i==0) {
				sMin[i] = 0.5 - trust;
				if(sMin[i] < infinityMin) sMin[i] = infinityMin;
				sMax[i] = 0.5 + trust;
				if(sMax[i]> infinityMax) sMax[i] = infinityMax;
			}
			else {
				sMin[i] = genome[i-1];      //range for design variables is fixed
				sMax[i] = genome[i-1];
			}
		}
		cout<<"after setting trust "<<endl;

		//get data points here
		
		if(currCount>neighCount) evalEuclidean = 1;
		else neighCount=currCount;
		if(evalEuclidean==1) {
   			for(i=0; i<currCount;i++) {
    				diffsq = 0;
				for(j=0; j<dimension;j++) {
					diffsq = diffsq +( ((inPoints[i][j+1]-genome[j])/(alphaMax[j]-alphaMin[j]))
								*((inPoints[i][j+1]-genome[j])/(alphaMax[j]-alphaMin[j])));
				}
				diffsq = sqrt(diffsq);
				inPoints[i][dimension+2] = diffsq;

			}
			cout<<"after get diffsq"<<endl;

			double dist[MAX_Indiv];
			double dist_index[MAX_Indiv];
			for(i=0; i<currCount; i++) {
				dist[i] = inPoints[i][dimension+2];
				dist_index[i] = i;
			}

			cout<<"before sorting"<<endl;
			sort(dist, dist_index, currCount, ASCENDING);
			cout<<"after sorting"<<endl;

			for(i=0; i<neighCount; i++) {
				index = (int)dist_index[i];
				for(j=0; j<dimension+2; j++) {
					dataPoints[i][j] = inPoints[index][j];
				}

			}
		}
		else {               // no need to measure the distance
			for(i=0; i<neighCount; i++) {
				for(j=0; j<dimension+2; j++) {
					dataPoints[i][j] = inPoints[i][j];
				}

			}
		}

		cout<<"after getting datapoint"<<endl;

		
		if(k==0) machInitial = 0.5;
		else {   //not first iteration
			if(rho>=0) {              //move towards right direction, so x(k) = x(k-1)
				machInitial = machApprox;
				fx = fxa; fax = faa;
			}
			else {fax = fx; }        // else wrong direction, still using the same point
		}
		
		
		float gen[25];
		gen[0] = machInitial;
		for(i=1; i<=24; i++) {
			gen[i] = genome[i-1];
		}

		cout<<"iter "<<k<<"before rbf"<<endl;
		//dimension + 1 ,   1 for mach
		RBFLocalSearch(neighCount, dimension+1, kerneltype, kernelhyp, Regparam, nf, ineqnf, eqnf, maxiter, sMin, sMax, dataPoints, gen, maximize);
		cout<<"iter "<<k<<"after rbf"<<endl;

		//get the result here , store it in these 2 vars below
		
		ifstream surrogateFile("surrogateOP.out");
		char dummy[40];
		for(i=0; i<9; i++) {
			surrogateFile>>dummy;
		}
		for(i=0; i<dimension+1; i++) {
			if(i==0) surrogateFile>>machApprox;
			else surrogateFile>>genomeApprox[i-1];
		}
		surrogateFile>>resultApprox;
		faa = resultApprox;

		//get fxa
		int j;
		char argu3[40]="./lana ";
		char argu4[20];
		sprintf(argu4,"%f 2.0 ", machApprox );
		strcat(argu3,argu4);
		ftyp.open("ex-alpha24.dat",ios::out|ios::trunc);
		for( j=0;j<dimension;j++) {
			ftyp<<genomeApprox[j]<<'\n';
		}
		ftyp.close();
		system(argu3);					// Calculate the exact fitness
		getPd.open("result.dat",ios::in);		// retrieve fitness output
		getPd>>pdOut;
		getPd.close();
		float fxa = pdOut;


		//store in database
		int repeatPoint = checkRepeatDesign(machApprox, genomeApprox);
  		if(!repeatPoint) {
  			currCountHash=currCount%MAX_Indiv;
  			inPoints[currCountHash][0] = machApprox;
  			for(j=0; j<dimension; j++) {
  				inPoints[currCountHash][j+1] = genomeApprox[j];
  			}
  			inPoints[currCountHash][dimension+1] = pdOut;
  			currCount++;
  		}

		
		if (k==1) fax = fx;  // for first iteration, approximation == exact
		
		
		
		rho = (fx - fxa) / (fax - faa);
		if(rho<=0.25) trustChange = 0.25;
		else if(rho>0.25 && rho<0.75) trustChange=1.0;
		else {
			float diff = machApprox - machInitial;
			if(diff<0.0) diff*=-1.0;
			if(diff==trust) trustChange = 2.0;
			else trustChange = 1.0;
		}
		//if(rho>=0) { fx = fxa; fax = faa; }
		//else { fax = fx; }

		trust = trustChange*trust;



	}

	return fxa;
}


int checkRepeatDesign(float machApprox, float* genomeApprox) {
	int i,j, repeat;	
	float difference;
	for(j=0;j<currCount;j++) {
		repeat=0;
		if (inPoints[j][0]==machApprox);
			repeat++;
		
		for( i=0;i<dimension;i++) {

			difference = inPoints[j][i+1] - genomeApprox[i];
			if (difference < 0)
				difference = -1*difference;
				
			if (difference < 0.0009)
				repeat+=1;

		}
	}
	if (repeat==dimension+1) return 1;
	else return 0;
}

////////////////////////////////////////////////////////////////////////////////////////////////


/* Functions Used by RBF Approximation */
void RBFLocalSearch(int traindatasize, int no_of_var, int kerneltype, float kernelhyp, float Regparam, int nf, int ineqnf, int eqnf, int maxiter, float *sminbounds, float *smaxbounds, float (*traindata)[dimension+2], float *startdesign, int maximize)
{
int w,k;
FILE *approxDatafp;

		approxDatafp = fopen("datafile.dat", "w");
		fprintf(approxDatafp,"%d  ! No. of training points\n",traindatasize);	// 500      ! No. of training points
		fprintf(approxDatafp,"%d  ! No. of input variables\n",no_of_var); 		// 2        ! No. of input variables
		/*
		ikernel ! Kernel type
          	0 -> Gaussian Kernel
          	1 -> Linear splines
          	2 -> Thin plate splines
          	3 -> Cubic splines
          	4 -> Multiquadrics
          	*/
		fprintf(approxDatafp,"%d  ! Kernel type\n",kerneltype); 				// 1        ! Kernel type
		fprintf(approxDatafp,"%f  ! kernel hyperparameter\n", kernelhyp); 		// 1.0      ! kernel hyperparameter
		fprintf(approxDatafp,"%f  ! Regularization parameter\n",Regparam); 		// 10E-10    ! Regularization parameter
		fprintf(approxDatafp,"%d  ! No. of objective functions\n",nf); 		// 1        ! No. of objective functions

		fprintf(approxDatafp,"%d  ! No. of inequality constraints\n",ineqnf); 	// 0        ! No. of inequality constraints
		fprintf(approxDatafp,"%d  ! No. of equality constraints\n",eqnf); 	// 0        ! No. of equality constraints

		fprintf(approxDatafp,"%d  ! print flag\n",0); 				// 0        ! print flag
		fprintf(approxDatafp,"%d  ! Maximum number of local optimization iterations\n",maxiter);	// 50       ! Maximum number of optimization iterations
		fprintf(approxDatafp,"%d  ! Mode of FFSQP operation\n",100); 		// 100      ! Mode of FFSQP operation

		for (w = 0; w < no_of_var; w++)
		{
			fprintf(approxDatafp,"%f  ",sminbounds[w]); 		// 0.0  0.0 ! Lower bounds
			cout<<"In RBFLocalSearch sminbounds["<<w<<"]:= "<<sminbounds[w]<<endl;
		}
		fprintf(approxDatafp,"  ! Lower bounds\n");

		for (w = 0; w < no_of_var; w++)
		{
			fprintf(approxDatafp,"%f  ",smaxbounds[w]); 		// 5.0  5.0 ! Upper bounds
			cout<<"In RBFLocalSearch smaxbounds["<<w<<"]:= "<<smaxbounds[w]<<endl;
		}
		fprintf(approxDatafp,"  ! Upper bounds\n");

		printf("Initial Guess!!!\n");
		for (w = 0; w < no_of_var; w++)
		{
			// for debugging
			//printf("%f  ",startdesign[w]);		 	// 2.2  2.9 ! Initial Guess
			fprintf(approxDatafp,"%f  ",startdesign[w]); 	// 2.2  2.9 ! Initial Guess
		}

		//printf("\n");
		fprintf(approxDatafp,"  ! Initial Guess\n");

		for (k = 0; k < traindatasize; k++)
		{
			for (w = 0; w < no_of_var; w++)
			{
				//printf("%f  ",traindata[k][w]);  // k+1 indicates the first nearest neighbour which is the same as the inital guess point will not be used for modelling
				fprintf(approxDatafp,"%e  ",traindata[k][w]);  // k+1 indicates the first nearest neighbour which is the same as the inital guess point will not be used for modelling
			}

			if (maximize == 0)
				fprintf(approxDatafp,"%e  ", traindata[k][no_of_var]); // write objective function value to datafile
			else
				fprintf(approxDatafp,"%e  ", -1*traindata[k][no_of_var]); // write objective function value to datafile

			for (w = 0; w < ineqnf+eqnf; w++)
			{
				fprintf(approxDatafp,"%e  ", traindata[k][no_of_var+w]); // write inequality and equality function values to datafile
			}

			fprintf(approxDatafp,"\n");
		}
		fclose(approxDatafp);
		//printf("Closing /staff/asysong/metaConfoilPressLoc/Model/datafile.dat!!!\n");
		//printf("Going for surrogate system call!!!\n");

		cout<<"=============================================================================="<<endl;
		system ("../pbn/locsearch < datafile.dat > surrogateOP.out");
		cout<<"=============================================================================="<<endl;
}
