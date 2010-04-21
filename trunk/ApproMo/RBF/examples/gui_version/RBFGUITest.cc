#include<iostream.h>
#include<time.h>
#include<fstream.h>
#include<vector.h>
#include "rbf.h"
#include "LinAlg/linalg.h"
#include "ArrayOp.h"
#define ELM 3

double sphere(vector<double> v) {
	double sum = 0;
	int i;
	int inDim = v.size();
	for(i=0; i<inDim;i++) {
		sum+=v[i]*v[i];
	}
	return sum;
}

double rosenbrock(vector<double> l_array)
{
	int varsize = l_array.size();
	int w;
	double l_value;
	
	l_value = 0;
	/* Rosenbrock's Function,Minimize f(x) = sum(100*(x(i)-x(i-1)^2)^2 + (1-x(i-1))^2) */
	// -3.0 < Bounds < 3.0
	for (w = 1; w < varsize; w++) {
		l_value = l_value + (100*((l_array[w] - (l_array[w-1]*l_array[w-1])) * (l_array[w] - (l_array[w-1]*l_array[w-1]))) + ((1 - l_array[w-1])*(1 - l_array[w-1])));
	}

	return l_value;
}

double rastrigin(vector<double> l_array)
{
	int varsize = l_array.size();
	int w,i,j,a=2;float c;
	double l_value, unNormVar;
  
  	
	l_value = 0;
	for (w = 0; w < varsize; w++) {
		l_value = l_value + (((l_array[w])*(l_array[w])) - 10.0*cos(2*M_PI*(l_array[w])));
	}
	l_value = (l_value + 10.0*varsize);
	
	return l_value;
}

double griewank(  vector< double > l_array )
{
        double l_value, l_Sumobj, l_Productobj;
        int w;
        int varsize = l_array.size();

        l_value = 0;
        l_Sumobj = 0;
        l_Productobj = 1;

        // Griewank test function, f11 in literature,Minimize f(x) = 1/4000*sum(x(i)-100)^2 - prod((x(i)-100)/sqrt(i)) + 1
        // -600.0 < Bounds < 600.0
        for (w = 0; w < varsize; w++) {
                //l_Sumobj = l_Sumobj + (pow(l_array[w],2)/4000);
                l_Sumobj = l_Sumobj + ((l_array[w] * l_array[w])/4000);
                l_Productobj = l_Productobj * cos(l_array[w]/(4*sqrt(w+1)));
        }
        l_value = (l_Sumobj + 1 - l_Productobj);
	return l_value;
}


double ackley(vector<double> l_array)
{
	int w;
	
	int varsize = l_array.size();
	double l_value, sqterm, costerm, a, b;
	
	l_value = 0;
	sqterm = 0;
	costerm = 0;
	
	/* Ackley's Function,f10 in literature f(x) = -20exp[-0.2*sqrt(1/n * sum(x(i)^2)) - exp(1/n * sum(cos2*pi*x(i)) + 20 + e */
	// -32.0 < Bounds < 32.0
	for (w = 0; w < varsize; w++) {
		sqterm = sqterm + (l_array[w]*l_array[w]);
		costerm = costerm + cos(2*M_PI*l_array[w]);
	}

	//l_value = -20*a - b + 20 + exp(1);
	l_value = ((-20)*exp((-0.2)*sqrt(sqterm/varsize))) - (exp(costerm/varsize)) + 20 + exp(1);
	//printf("l_value:=%f\n",l_value);
	return l_value;
}


int main(int argc, char** argv) {
	int i,j,k;
	//RBF aRBF( 5, 1, gaussian, 0.4, 0.5);
	//RBF aRBF();

	char lm[50]; learnMethod lmU;
	char tf[20];
	unsigned inDim = 2,  nInput;
	char kernFunc[40]; kernelType kernFuncU;
	char errorMeth[5]; errMethod errorMethU;
	double loBound=-4, wholeRange = 8, kernWidth, maxError, regParam;
	unsigned maxBasis, nBasis, maxClusterIter, maxRegOptIter;

	char activation[20]; unsigned activationU;
	unsigned maxNHidden;

	double (*evaluate)(vector<double>);
	
	sprintf(lm,"%s",argv[1]);
	if(!strcmp(lm,"Interpolate")) { lmU=interpolate; }
	else if(!strcmp(lm,"Regression")) { lmU=ridgeRegression; }
	else if(!strcmp(lm,"OLSForward")) { lmU=OLS_ForwardSelection; }
	else {lmU = ELM;}

	sprintf(tf,"%s",argv[2]);
	if(!strcmp(tf,"ackley")) {evaluate = &ackley;}
	else if(!strcmp(tf,"griewank")) {evaluate = &griewank;}
	else if(!strcmp(tf,"rosenbrock")) {evaluate = &rosenbrock;}	
	else if(!strcmp(tf,"sphere")) {evaluate = &sphere;}
	else if(!strcmp(tf,"rastrigin")) {evaluate = &rastrigin;}

	nInput = atoi(argv[3]);
	
	//rbf learning methods
	if(lmU<ELM) {
		

		sprintf(kernFunc,"%s",argv[4]);
		if(!strcmp(kernFunc,"gaussian")) {kernFuncU = gaussian;}
		else if(!strcmp(kernFunc,"linear_spline")) {kernFuncU = linear_spline;}
		else if(!strcmp(kernFunc,"cubic_spline")) {kernFuncU = cubic_spline;}
		else if(!strcmp(kernFunc,"thin_plate_spline")) {kernFuncU = thin_plate_spline;}
		else if(!strcmp(kernFunc,"multiquadric")) {kernFuncU = multiquadric;}
		else if(!strcmp(kernFunc,"inverse_multiquadric")) {kernFuncU = inverse_multiquadric;}
		else if(!strcmp(kernFunc,"cauchy")) {kernFuncU = cauchy;}

		kernWidth = atof(argv[5]);
		
		if(lmU==ridgeRegression) {
			nBasis = atoi(argv[6]);
			regParam = atof(argv[7]);
			maxError = atof(argv[8]);
			sprintf(errorMeth,"%s",argv[9]);

			if(!strcmp(errorMeth,"GCV")) {errorMethU=GCV;}
			else if(!strcmp(errorMeth,"FPE")) {errorMethU=FPE;}
			else if(!strcmp(errorMeth,"UEV")) {errorMethU=UEV;}
			else if(!strcmp(errorMeth,"BIC")) {errorMethU=BIC;}

			maxClusterIter = atoi(argv[10]);
			maxRegOptIter = atoi(argv[11]);
		}
		else if(lmU==OLS_ForwardSelection) {
			
			maxBasis = atoi(argv[6]);
			maxError = atof(argv[7]);
			
		}
	}
	else {		//it is ELM
		sprintf(activation,"%s",argv[4]);
		maxNHidden = atoi(argv[5]);
		maxError = atof(argv[6]);
	}

	

	Array<double> inArray(nInput, inDim);
	Array<double> targArray(nInput, 1);
	
	//might need to change to lhs later
	srand(time(NULL));
	for(i=0; i<nInput; i++) {
		for(j=0; j<inDim; j++) 
			inArray(i,j) = loBound+wholeRange*(rand()%12134)/12133.0;
	}
		
	for(i=0; i<nInput; i++) {
		vector<double> arr(inDim);
		for(j=0; j<inDim; j++) {
			arr[j] = inArray(i,j);
		}
		targArray(i,0) = evaluate(arr);
	}
// 	for(i=0; i<nInput; i++) {
// 		for(j=0; j<inDim; j++) 
// 			cout<<inArray(i,j)<<"\t";
// 		//cout<<targArray(i,0)<<"\t"<<endl;
// 		cout<<endl;
// 	}

	RBF aRBF;
	if(lmU<ELM) {
		Array<double> kernParams(1,1);
		kernParams(0,0) = kernWidth;
	
		
		switch(lmU) {
			case interpolate: aRBF = RBF(inDim,  kernFuncU, kernParams);
				break;
			case ridgeRegression: aRBF = RBF(inDim, kernFuncU, kernParams, nBasis, regParam, maxError, errorMethU, maxClusterIter, maxRegOptIter);
				break;
			case OLS_ForwardSelection: aRBF = RBF(inDim, kernFuncU, kernParams, maxBasis,  maxError);
				break;
			default: break;
		}
		
	}

	aRBF.train(inArray, targArray);

	return 0;
}


