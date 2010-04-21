#include "rbf.h"
#include <iostream.h>
#include "sqr.h"
#include "LinAlg/linalg.h"
#include <fstream.h>
#include <iomanip.h>



RBF::RBF() {
  
}

RBF::~RBF() {
  
}


RBF::RBF( unsigned inDim, kernelType aKernel, Array<double> kernParams) {   //interpolation
	setParams(inDim, aKernel, kernParams);
	learner = interpolate;
	regularParam = 0;
}

RBF::RBF( unsigned inDim, kernelType aKernel, Array<double> kernParams, unsigned nKernel, double regParam, double maxErr, errMethod anErrMethod, unsigned maxKMeanIter, unsigned maxRegIter) {	//supervised + ridge regression
	setParams(inDim, aKernel, kernParams);
	learner = ridgeRegression;
	regularParam = regParam;
	nBasis = nKernel;
	maxError = maxErr;
	errorMethod = anErrMethod;
	maxKMeanIteration = maxKMeanIter;
	maxRegIteration = maxRegIter;

	switch(anErrMethod) {
		case GCV: errFuncPtr = &GCVFunc; estRegParamPtr = &estRegParamGCV;
			break;
		case UEV: errFuncPtr = &UEVFunc; estRegParamPtr = &estRegParamUEV;
			break;		
		case FPE: errFuncPtr = &FPEFunc; estRegParamPtr = &estRegParamFPE;
			break;		
		case BIC: errFuncPtr = &BICFunc; estRegParamPtr = &estRegParamBIC;
			break;
		default: errFuncPtr = &GCVFunc; estRegParamPtr = &estRegParamGCV;
			break;
	}
}

RBF::RBF( unsigned inDim, kernelType aKernel, Array<double> kernParams, unsigned maximumBasis, double maxErr) {	//unsupervised , OLS-forward selection
	setParams(inDim, aKernel, kernParams);
	maxError = maxErr;
	maxNBasis = maximumBasis;
	learner=OLS_ForwardSelection;
}

void RBF::setParams(unsigned inDim, kernelType aKernel, Array<double> kernParams) {
	inputDim = inDim;
	kernelFunction = aKernel;
	kernelParams = kernParams;
	
	switch(kernelFunction) {
		case gaussian: kernelFuncPtr = &gaussFunc;
			break;
		case linear_spline: kernelFuncPtr = &linearFunc;
			break;		
		case cubic_spline: kernelFuncPtr = &cubicFunc;
			break;		
		case thin_plate_spline: kernelFuncPtr = &thinPlateFunc;
			break;
		case multiquadric: kernelFuncPtr = &multiquadricFunc;
			break;
		case inverse_multiquadric: kernelFuncPtr = &inverseMultiquadricFunc;
			break;
		case cauchy: kernelFuncPtr = &cauchyFunc;
			break;
		default: kernelFuncPtr = &gaussFunc;
			break;
	}
}

//=======================================================================

double RBF::mse( Database data ) {

	Array < double > inputData(1,1);
	Array < double > targetData(1,1);
	data.GetLastData(inputData,targetData);
	
	Array < double > outputData(1,1);
	evaluate ( inputData, outputData );
	
	double rms_LB = 0;
	for ( int i = 0; i < (int)inputData.dim(0); i++) {
	rms_LB += sqr( outputData(i,0) - targetData(i,0) );
	}
	
	return (rms_LB/(double)inputData.dim(0));
}

//=======================================================================

double RBF::mse( Array < double > input,
		     Array < double > target ) {
	
	Array < double > inputData(1,1);
	Array < double > targetData(1,1);
	inputData.resize ( input.dim (0), input.dim(1) );
	targetData.resize ( target.dim(0), target.dim(1) );
	for ( int i = 0; i < (int)input.dim(0); i++ ) {
	for ( int j = 0; j < (int)input.dim(1); j++ ) {
	inputData (i,j) = input (i,j);
	
	}
		targetData (i,0) = target (i,0);
	}
	
	
	Array < double > outputData(inputData.dim(0),1);
	
		
	evaluate ( inputData, outputData );
	
		
	double rms_LB = 0;
	for ( int i = 0; i < (int)inputData.dim(0); i++) {
	rms_LB += sqr( outputData(i,0) - targetData(i,0) );
	}
	
	return (rms_LB/(double)inputData.dim(0));
}

//=======================================================================

double RBF::mse( Population offsprings ) {

	// the parameters have to be in the first Chromosome!!!
	unsigned inputDimension = dynamic_cast < ChromosomeT < double >& > (offsprings[0][0]).size();
	unsigned numberOfOffsprings = offsprings.size();
	
	Array < double > inputData (numberOfOffsprings,inputDimension);
	Array < double > targetData (numberOfOffsprings, 1);
	
	for ( unsigned i = 0; i < numberOfOffsprings; i++ ) {
	for ( unsigned j = 0; j < inputDimension; j++ ) {
	inputData(i,j) = dynamic_cast < ChromosomeT < double >& > (offsprings[i][0] )[j];
	}
	targetData (i,0) = offsprings[i].fitnessValue();
	}  
	
	double error = mse ( inputData, targetData );
	
	return error;
}



int RBF::train( Array < double > inputData, Array < double > targetData ) {
	unsigned i,j,k;
	if(learner==interpolate) {
		nBasis = inputData.dim(0);
		weightMatrix.resize(nBasis, 1, true);
		designMatrix.resize(inputData.dim(0), nBasis, true);
		baseCentres.resize(nBasis, inputDim, true);
		initArrayToZero(weightMatrix);
		initArrayToZero(designMatrix);
		initArrayToZero(baseCentres);
	
		for(i=0; i<nBasis; i++) {
			for(j=0; j<inputDim; j++) {
				baseCentres(i,j) = inputData(i,j);		
			}
		}
		ofstream outFile("./plotInterpolate/outInterpolate.dat",ios::trunc);
		outFile<<"Setting design matrix..."<<endl;
		cout<<"Setting design matrix..."<<endl;
		outFile.close();
			
		setDesignMatrix(inputData);
		
		outFile.open("./plotInterpolate/outInterpolate.dat",ios::app);
		outFile<<"Setting weight matrix..."<<endl;
		cout<<"Setting weight matrix..."<<endl;
		outFile.close();

		setWeightMatrix(targetData);
		
		outFile.open("./plotInterpolate/outInterpolate.dat",ios::app);
		outFile<<"Calculating test points..."<<endl;
		cout<<"Calculating test points..."<<endl;
		outFile.close();	
	
		printToFile(0, weightMatrix, baseCentres, kernelParams, kernelFunction, learner);   //remove this line in
		
		outFile.open("./plotInterpolate/outInterpolate.dat",ios::app);
		outFile<<"Done..."<<endl;
		cout<<"Done..."<<endl;
		outFile.close();
	}
	else if(learner==ridgeRegression) {
		weightMatrix.resize(nBasis, 1, true);
		designMatrix.resize(inputData.dim(0), nBasis, true);
		baseCentres.resize(nBasis, inputDim, true);
		initArrayToZero(weightMatrix);
		initArrayToZero(designMatrix);
		initArrayToZero(baseCentres);

		ofstream outFile("./plotRegression/outRegression.dat",ios::trunc);
		if(nBasis==inputData.dim(0)) {
			cout<<"K-Mean algorithm is not necessary for the interpolation case. Proceed with learning..."<<endl;
			outFile<<"K-Mean algorithm is not necessary for the interpolation case. Proceed with learning..."<<endl;
			for(i=0; i<nBasis; i++) {
				for(j=0; j<inputDim; j++) {
					baseCentres(i,j) = inputData(i,j);		
				}
			}
		}
		else {
			Array<double> centresFound(nBasis, inputDim);
			Array<unsigned> clustMap(inputDim,1);
			KMean(nBasis, maxKMeanIteration, inputData, centresFound, clustMap, outFile);
			for(i=0; i<nBasis; i++) {
				for(j=0; j<inputDim; j++) {
					baseCentres(i,j) = centresFound(i,j);		
				}
			}
		}
		outFile.close();
	
		setDesignMatrix(inputData);
		setWeightMatrix(targetData);
		double error=HUGE;

		char errString[5];
		switch (errorMethod) {
			case GCV: sprintf(errString, "GCV"); break;
			case FPE: sprintf(errString, "FPE"); break;
			case UEV: sprintf(errString, "UEV"); break;
			case BIC: sprintf(errString, "BIC"); break;
		}
		
		cout<<endl<<"Start ridge regression..."<<endl;
		cout<<"iteration	regParam	        "<<errString<<" error	        mean square error"<<endl;
		cout<<"-----------------------------------------------------------------------------------------"<<endl;

		outFile.open("./plotRegression/outRegression.dat",ios::app);
		outFile<<endl<<"Start ridge regression..."<<endl;
		outFile<<"iteration	regParam	        "<<errString<<" error	        mean square error"<<endl;
		outFile<<"---------------------------------------------------------------------------------------"<<endl;
		
		for(i=0; i<maxRegIteration && error>maxError; i++) {
			setWeightMatrix(targetData);
			Array<double> output(inputData.dim(0),1);
			evaluate(inputData, output);

			outFile.close();
			printToFile(i, weightMatrix, baseCentres, kernelParams, kernelFunction, learner);   //remove this line in later version
			
			double optError = errFuncPtr(designMatrix, output, regularParam);
			error=mse(inputData, targetData);	

			cout<<std::right<<i<<"\t\t"<<setw(10)<<std::left<<regularParam<<"\t\t"<<setw(10)<<std::left<<optError<<"\t\t"<<setw(10)<<std::left<<error<<endl;
			
			outFile.open("./plotRegression/outRegression.dat",ios::app);
			outFile<<std::right<<i<<"\t\t"<<setw(10)<<std::left<<regularParam<<"\t\t"<<setw(10)<<std::left<<optError<<"\t\t"<<setw(10)<<std::left<<error<<endl;
			
			regularParam = estRegParamPtr(regularParam, weightMatrix, designMatrix, output);
		}
		outFile<<"Done..."<<endl;
		cout<<"Done..."<<endl;
		outFile.close();
 	}
	else if(learner==OLS_ForwardSelection) {

		unsigned nInput = inputData.dim(0);
		Array<unsigned> selectedMap(nInput,1);
		initArrayToZero(selectedMap);

		
		nBasis=1;	//current number of basis function
		double error=HUGE;
		
		
		orthogonalMatrix.resize(nInput,nBasis,true);
		initArrayToZero(orthogonalMatrix);

		orthogonalWeightMatrix.resize(nBasis,1,true);
		initArrayToZero(orthogonalWeightMatrix);

		designMatrix.resize(nInput,nBasis,true);
		initArrayToZero(designMatrix);

		// iteration 1, try all input to get the first centre
		unsigned maxIndex;
		double errReduction;
		double maxErrReduction = -1;
		Array<double> newOrtWeight(1,1);
		Array<double> tempOrtWeight(1,1);
		Array<double> newOrtVector(nInput, 1);
		Array<double> tempOrtVector(nInput, 1);
		Array<double> newA(nBasis, 1);
		Array<double> tempA(nBasis, 1);
	
		Array<double> temp, temp2;
		Array<double> tempBaseCentre(1,inputDim);
		Array<double> newBaseCentre(1,inputDim);


		ofstream outFile("./plotOLSForward/outOLSForward.dat",ios::trunc);
		cout<<endl<<"Start Orthogonal Least Square Forward Selection..."<<endl<<endl;
		cout<<"iteration	node added		mean square error"<<endl;
		cout<<"-----------------------------------------------------------------------------------------"<<endl;

		outFile<<endl<<"Start Orthogonal Least Square Forward Selection..."<<endl<<endl;
		outFile<<"iteration	node added		 mean square error"<<endl;
		outFile<<"--------------------------------------------------------------------------------------"<<endl;

		//cout<<"try to get first node"<<endl;
		//outFile<<"try to get first node"<<endl;
		for(i=0; i<nInput; i++) {
			//cout<<"try node "<<i<<endl;

			for(j=0; j<inputDim; j++) tempBaseCentre(0,j) = inputData(i,j);
			baseCentres = tempBaseCentre;

			setDesignVector(inputData, 0, baseCentres, tempOrtVector);	
		
		
			temp = innerProduct( transpose(tempOrtVector), targetData );
			temp2 = innerProduct( transpose(tempOrtVector), tempOrtVector );

		
			tempOrtWeight(0,0) = temp(0,0)/temp2(0,0);
	
			temp = innerProduct( transpose(targetData), targetData );
			errReduction = sqr(tempOrtWeight(0,0))*temp2(0,0)/temp(0,0);
		
			if(i==0) {
				maxIndex = 0;
				maxErrReduction = errReduction;
				newOrtVector = tempOrtVector;
				newOrtWeight = tempOrtWeight;
				newBaseCentre = tempBaseCentre;
			//	cout<<"current best node is node no. "<<i<<endl;
			}
			else {
				if(errReduction>maxErrReduction) { 
					maxErrReduction = errReduction;
					maxIndex = i;
					newOrtVector = tempOrtVector;
					newOrtWeight = tempOrtWeight;	
					newBaseCentre = tempBaseCentre;
			//		cout<<"current best node is node no. "<<i<<endl;
				}
			}
		}
		
		//add the first chosen centre to design matrix H
		//cout<<"So, first node is node no. "<<maxIndex<<endl;
		designMatrix = newOrtVector;
		orthogonalMatrix = newOrtVector;
		orthogonalWeightMatrix = newOrtWeight;	
		baseCentres = newBaseCentre;
		selectedMap(maxIndex,0) = 1;	

		setWeightMatrix(targetData);
		error = mse(inputData, targetData);
		
	
		cout<<std::right<<0<<"\t\t"<<setw(10)<<std::left<<maxIndex<<"\t\t"<<setw(10)<<std::left<<error<<endl;
		outFile<<std::right<<0<<"\t\t"<<setw(10)<<std::left<<maxIndex<<"\t\t"<<setw(10)<<std::left<<error<<endl;

		printToFile(0, weightMatrix, baseCentres, kernelParams, kernelFunction, learner);   //remove this line in
		
		Array<double> tempDesignVector(nInput,1);
		Array<double> newDesignVector(nInput,1);

		unsigned t;
		//cout<<"try to get other nodes "<<endl;
		for(k=1; k<maxNBasis; k++) {
			if(error<maxError) { break; }
			
			nBasis++;

			tempA.resize(nBasis,1,true);

// 			cout<<"selected map :\n";
// 			for(j=0; j<nInput; j++) cout<<j<<"  "<<selectedMap(j,0)<<"\n";			

			bool first = true;
			for(i=0; i<nInput; i++) {	
				//find tempA
				if(selectedMap(i,0)==1) {  continue; }				

				Array<double> copyOfBaseCentres = baseCentres;

				for(j=0; j<inputDim; j++) tempBaseCentre(0,j) = inputData(i,j);
				copyOfBaseCentres = copyOfBaseCentres.append_rows(tempBaseCentre);

				
				setDesignVector(inputData, k, copyOfBaseCentres, tempDesignVector);
				for(j=0; j<nBasis-1; j++) {
					Array<double> qj(nInput,1);
					temp = transpose(orthogonalMatrix);
					for(t=0; t<nInput; t++) qj(t,0) = temp(j,t);
					temp = innerProduct( transpose(qj), tempDesignVector);
					temp2 = innerProduct( transpose(qj), qj);
					tempA(j,0) = temp(0,0)/temp2(0,0);
				}
				tempA(j,0) = 1;

				//find tempOrtVector
				temp.resize(nInput,1,false); initArrayToZero(temp);
				Array<double> orthoTransposed = transpose(orthogonalMatrix);
				for(t=0; t<nInput; t++) {
					for(j=0; j<nBasis-1; j++) {
						temp(t,0) += tempA(j,0)*orthoTransposed( j,t ); 
					}
				}				
				
				for(j=0; j<nInput; j++) {
					tempOrtVector(j,0) = tempDesignVector(j,0)-temp(j,0); 
				}

				//find max error reduction using tempOrtVector and targetData
				temp = innerProduct( transpose(tempOrtVector), targetData );
				temp2 = innerProduct( transpose(tempOrtVector), tempOrtVector );
	
				tempOrtWeight(0,0) = temp(0,0)/temp2(0,0);
		
				temp = innerProduct( transpose(targetData), targetData );
				errReduction = sqr(tempOrtWeight(0,0))*temp2(0,0)/temp(0,0);
			
				if(first) {
					maxIndex = i;
					maxErrReduction = errReduction;
					newDesignVector = tempDesignVector;
					newOrtVector = tempOrtVector;
					newOrtWeight = tempOrtWeight;
					newBaseCentre = tempBaseCentre;
					//newA =tempA;
				}
				else if(errReduction>maxErrReduction) { 
					maxErrReduction = errReduction;
					newDesignVector = tempDesignVector;
					maxIndex = i;
					newOrtVector = tempOrtVector;
					newOrtWeight = tempOrtWeight;	
					newBaseCentre = tempBaseCentre;	
					//newA = tempA;
				}	
				first = false;
			}

		//	cout<<"next node ( "<<k<<" ) is node no. "<<maxIndex<<endl; 
			//A = A.append_cols(newA);
			baseCentres = baseCentres.append_rows(newBaseCentre);
			designMatrix = designMatrix.append_cols(newDesignVector);
			orthogonalMatrix = orthogonalMatrix.append_cols(newOrtVector);	
			orthogonalWeightMatrix = orthogonalWeightMatrix.append_rows(newOrtWeight);
			selectedMap(maxIndex,0) = 1;

			setWeightMatrix(targetData);
			error = mse(inputData, targetData);

	
			cout<<std::right<<k<<"\t\t"<<setw(10)<<std::left<<maxIndex<<"\t\t"<<setw(10)<<std::left<<error<<endl;
			outFile<<std::right<<k<<"\t\t"<<setw(10)<<std::left<<maxIndex<<"\t\t"<<setw(10)<<std::left<<error<<endl;

			printToFile(k, weightMatrix, baseCentres, kernelParams, kernelFunction, learner);   //remove this line in
		} 

		outFile<<"Done..."<<endl;
		cout<<"Done..."<<endl;
		outFile.close();
		

	}

	ofstream finishFile;
	switch(learner) {
		case interpolate: finishFile.open("./plotInterpolate/finish.dat", ios::trunc); finishFile<<"0"<<endl;
			break;
		case ridgeRegression: finishFile.open("./plotRegression/finish.dat", ios::trunc); finishFile<<i-1<<endl; 
			break;
		case OLS_ForwardSelection: finishFile.open("./plotOLSForward/finish.dat", ios::trunc); finishFile<<k-1<<endl;
			break;
	}
	finishFile.close();

	return 0;
}

static void RBF::printToFile(unsigned iter, Array<double> wMatrix, Array<double> centres, Array<double> kernParams, kernelType kernFunction, learnMethod aMethod) {
	unsigned i,j, npts = 20;  //npts per dimension actually
	double wholeRange = 8.0;
	double loBound = -4;
	Array<double> grid(npts*npts,2);
	Array<double> gridOutput(npts*npts,1);
	for(i=0; i<npts; i++) {
		double a=i*wholeRange/npts+loBound;
		for(j=0; j<npts; j++) {
			double b=j*wholeRange/npts+loBound;
			grid(i*npts+j,0) = a;
			grid(i*npts+j,1) = b;
		}
	}
	char fileName[100];
	switch(aMethod) {
		case interpolate: sprintf(fileName, "./plotInterpolate/plot-%d.dat", iter);
			break;
		case ridgeRegression: sprintf(fileName, "./plotRegression/plot-%d.dat", iter);
			break;
		case OLS_ForwardSelection: sprintf(fileName, "./plotOLSForward/plot-%d.dat", iter);
			break;
	}

	ofstream outFile(fileName,ios::trunc);

	double(*localFuncPtr)(Array<double>, unsigned, Array<double>, unsigned, Array<double>);
	switch(kernFunction) {
		case gaussian: localFuncPtr = &gaussFunc;
			break;
		case linear_spline: localFuncPtr = &linearFunc;
			break;
		case cubic_spline: localFuncPtr = &linearFunc;
			break;
		case thin_plate_spline: localFuncPtr = &thinPlateFunc;
			break;
		case multiquadric: localFuncPtr = &multiquadricFunc;
			break;
		case inverse_multiquadric: localFuncPtr = &inverseMultiquadricFunc;
			break;
		case cauchy: localFuncPtr = &cauchyFunc;
			break;
		default: localFuncPtr = &gaussFunc;
			break;
	}


	Array<double> desMatrix(grid.dim(0), centres.dim(0));
	for(i=0; i<desMatrix.dim(0); i++) {
		for(j=0; j<desMatrix.dim(1); j++) {
			desMatrix(i,j) = localFuncPtr(grid, i, centres, j, kernParams);
		}
	}
	
	gridOutput = innerProduct(desMatrix, wMatrix);



	for(i=0; i<npts; i++) {
		for(j=0; j<npts; j++) {
			outFile<<grid(i*npts+j,0)<<"\t"<<grid(i*npts+j,1)<<"\t"<<gridOutput(i*npts+j,0)<<endl;
		}
		outFile<<endl;
	}
	
	char fileName2[100];
	switch(aMethod) {
		case interpolate: sprintf(fileName2, "./plotInterpolate/plot-%d.tcl", iter);
			break;
		case ridgeRegression: sprintf(fileName2, "./plotRegression/plot-%d.tcl", iter);
			break;
		case OLS_ForwardSelection: sprintf(fileName2, "./plotOLSForward/plot-%d.tcl", iter);
			break;
	}
	
	ofstream gnuplotFile("gnuplot.dat", ios::trunc);
	gnuplotFile<<"set term tkcanvas"<<endl
		<<"set output \'"<<fileName2<<"\'"<<endl
		<<"splot \'"<<fileName<<"\' w line"<<endl;
	
	system("/usr/bin/gnuplot 'gnuplot.dat'");
}

void RBF::setWeightMatrix(Array<double>& targetData) {
	Array<double> A = innerProduct(transpose(designMatrix),designMatrix);

	unsigned i,j;

	if(learner==ridgeRegression) {
		Array<double> identity(nBasis, nBasis);
		initArrayToZero(identity);
	
		for(i=0; i<nBasis; i++) {
			for(j=0; j<nBasis; j++) {
				if(i==j) identity(i,j) = regularParam;
				A(i,j)+=identity(i,j);
			}
		}
	}
	A = invert(A);
	
	weightMatrix = innerProduct(innerProduct(A, transpose(designMatrix)), targetData); 
}

void RBF::setDesignMatrix(Array<double>& inputData) {
	unsigned i, j;

	
	for(i=0; i<designMatrix.dim(0); i++) {
		for(j=0; j<designMatrix.dim(1); j++) {
			designMatrix(i,j) = kernelFuncPtr(inputData, i, baseCentres, j, kernelParams);
		}
	}

}

void RBF::setDesignVector(Array<double> input, unsigned index, Array<double> centres, Array<double>& outputVect) {
	unsigned i, j;
	

	for(i=0; i<input.dim(0); i++) {
		outputVect(i,0) = kernelFuncPtr(input, i, centres, index, kernelParams);
	}

}

static double RBF::gaussFunc(Array<double> input, unsigned i, Array<double> centre, unsigned j, Array<double> params) {
	//double sigma = 1.0;   //change later
	double sigmaParam = params(0,0);
	double distance = getDistance(input,i,centre,j);
	
	return exp(-1.0*(distance/sigmaParam)*(distance/sigmaParam));
}

static double RBF::cubicFunc(Array<double> input, unsigned i, Array<double> centre, unsigned j, Array<double> params) {
	double distance = getDistance(input, i, centre, j);
	return sqr(distance)*distance;
}

static double RBF::linearFunc(Array<double> input, unsigned i, Array<double> centre, unsigned j, Array<double> params) {
	
	return  getDistance(input, i, centre, j);
}

static double RBF::thinPlateFunc(Array<double> input, unsigned i, Array<double> centre, unsigned j, Array<double> params) {
	double distance = getDistance(input, i, centre, j);
	//get the code later
	return distance;
}

static double RBF::multiquadricFunc(Array<double> input, unsigned i, Array<double> centre, unsigned j, Array<double> params) {
	double distance = getDistance(input, i, centre, j);
	return sqrt(sqr(distance)+sqr(params(0,0)))/params(0,0);
}

static double RBF::inverseMultiquadricFunc(Array<double> input, unsigned i, Array<double> centre, unsigned j, Array<double> params) {
	double distance = getDistance(input, i, centre, j);
	return params(0,0)/sqrt(sqr(distance)+sqr(params(0,0)));
}

static double RBF::cauchyFunc(Array<double> input, unsigned i, Array<double> centre, unsigned j, Array<double> params) {
	double distance = getDistance(input, i, centre, j);
	return sqr(params(0,0))/(sqr(distance)+sqr(params(0,0)));
}



static double RBF::YPYFunc(Array<double> desMatrix, Array<double> output, double regParam, Array<double>& A) {

	Array<double> temp = innerProduct(transpose(desMatrix),desMatrix);
	Array<double> identity(desMatrix.dim(1),desMatrix.dim(1));



	unsigned i,j;

	for(i=0; i<identity.dim(0); i++) {
		for(j=0; j<identity.dim(1); j++) {
			if(i==j) identity(i,j) = regParam;
			else identity(i,j) = 0;
		}
	}
	
	for(i=0; i<identity.dim(0); i++) {
		for(j=0; j<identity.dim(1); j++) {
			temp(i,j)+=identity(i,j);		
		}
	}
		
	A = invert(temp);

	Array<double> PY = innerProduct(innerProduct(innerProduct(desMatrix,A), transpose(desMatrix)), output) ;
	for(i=0; i<PY.dim(0); i++) {
		for(j=0; j<PY.dim(1); j++) {
			PY(i,j) = output(i,j) - PY(i,j);
		}
	}

	return trace( innerProduct(transpose(PY), PY) );
}

static double RBF::GCVFunc(Array<double> desMatrix, Array<double> output, double regParam) { // Generalized Cross-Validation
	double m = desMatrix.dim(1);
	double p = output.dim(0);
	Array<double> A(m,m);
	double ypy = YPYFunc(desMatrix, output, regParam, A);
	double g = m-regParam*trace(A);
	
	double psi = sqr(p/(p-g));

	return psi*ypy/p;
}

static double RBF::UEVFunc(Array<double> desMatrix, Array<double> output, double regParam) { // Unbiased Estimate Variance
	double m = desMatrix.dim(1);
	double p = output.dim(0);
	Array<double> A(m,m);
	double ypy = YPYFunc(desMatrix, output, regParam, A);
	double g = m-regParam*trace(A);
	
	double psi = p/(p-g);

	return psi*ypy/p;}

static double RBF::FPEFunc(Array<double> desMatrix, Array<double> output, double regParam) { // Final Prediction Error
	double m = desMatrix.dim(1);
	double p = output.dim(0);
	Array<double> A(m,m);
	double ypy = YPYFunc(desMatrix, output, regParam, A);
	double g = m-regParam*trace(A);
	
	double psi = (p+g)/(p-g);

	return psi*ypy/p;
}

static double RBF::BICFunc(Array<double> desMatrix, Array<double> output, double regParam) { // Bayesian Information Criterion
	double m = desMatrix.dim(1);
	double p = output.dim(0);
	Array<double> A(m,m);

	
	double ypy = YPYFunc(desMatrix, output, regParam, A);
	
	
	double g = m-regParam*trace(A);
	
	double psi = (p+g*(log(p)-1))/(p-g);

	return psi*ypy/p;
}

static double RBF::estRegParamUEV(double regParam, Array<double> wMatrix, Array<double> desMatrix, Array<double>output) {
	unsigned i,j;

	double m = desMatrix.dim(1);
	double p = output.dim(0);
	Array<double> A(m,m);

	
	double ypy = YPYFunc(desMatrix, output, regParam, A);
	double g = m-regParam*trace(A);

	Array<double> A2 = innerProduct(A,A);
	Array<double> A3 = innerProduct(A2,A);

	Array<double> temp1(m,m);
	for(i=0; i<m; i++) {
		for(j=0; j<m; j++) {
			temp1(i,j) = A(i,j)-regParam*A2(i,j);	
		}
	}

	Array<double> HY = innerProduct(transpose(desMatrix),output);
	Array<double> temp2 = innerProduct( innerProduct( transpose(HY), A3 ), HY);
	

	double eta = 1.0/(2*(p-g));
	return eta * ypy * trace(temp1) / trace(temp2);
}

static double RBF::estRegParamGCV(double regParam, Array<double> wMatrix, Array<double> desMatrix, Array<double>output) {
		unsigned i,j;

	double m = desMatrix.dim(1);
	double p = output.dim(0);
	Array<double> A(m,m);

	
	double ypy = YPYFunc(desMatrix, output, regParam, A);
	double g = m-regParam*trace(A);

	Array<double> A2 = innerProduct(A,A);
	Array<double> A3 = innerProduct(A2,A);

	Array<double> temp1(m,m);
	for(i=0; i<m; i++) {
		for(j=0; j<m; j++) {
			temp1(i,j) = A(i,j)-regParam*A2(i,j);	
		}
	}

	Array<double> HY = innerProduct(transpose(desMatrix),output);
	Array<double> temp2 = innerProduct( innerProduct( transpose(HY), A3 ), HY);
	

	double eta = 1.0/(p-g);
	return eta * ypy * trace(temp1) / trace(temp2);
}

static double RBF::estRegParamFPE(double regParam, Array<double> wMatrix, Array<double> desMatrix, Array<double>output) {
	unsigned i,j;

	double m = desMatrix.dim(1);
	double p = output.dim(0);
	Array<double> A(m,m);

	
	double ypy = YPYFunc(desMatrix, output, regParam, A);
	double g = m-regParam*trace(A);

	Array<double> A2 = innerProduct(A,A);
	Array<double> A3 = innerProduct(A2,A);

	Array<double> temp1(m,m);
	for(i=0; i<m; i++) {
		for(j=0; j<m; j++) {
			temp1(i,j) = A(i,j)-regParam*A2(i,j);	
		}
	}

	Array<double> HY = innerProduct(transpose(desMatrix),output);
	Array<double> temp2 = innerProduct( innerProduct( transpose(HY), A3 ), HY);
	

	double eta = p/((p-g)*(p+g));
	return eta * ypy * trace(temp1) / trace(temp2);
}

static double RBF::estRegParamBIC(double regParam, Array<double> wMatrix, Array<double> desMatrix, Array<double>output) {
	unsigned i,j;

	double m = desMatrix.dim(1);
	double p = output.dim(0);
	Array<double> A(m,m);

	
	double ypy = YPYFunc(desMatrix, output, regParam, A);
	double g = m-regParam*trace(A);

	Array<double> A2 = innerProduct(A,A);
	Array<double> A3 = innerProduct(A2,A);

	Array<double> temp1(m,m);
	for(i=0; i<m; i++) {
		for(j=0; j<m; j++) {
			temp1(i,j) = A(i,j)-regParam*A2(i,j);	
		}
	}

	Array<double> HY = innerProduct(transpose(desMatrix),output);
	Array<double> temp2 = innerProduct( innerProduct( transpose(HY), A3 ), HY);
	

	double eta = p*log(p)/( 2 * ( p - g ) * ( p + ( log(p) - 1 ) * g ) );
	return eta * ypy * trace(temp1) / trace(temp2);
}

Array<double> RBF::getDesignMatrix() {
	return designMatrix;
}

Array<double> RBF::getWeightMatrix() {
	return weightMatrix;
}

Array<double> RBF::getBaseCentres() {
	return baseCentres;
}

Array<double> RBF::getOrthogonalMatrix() {
	if(learner==OLS_ForwardSelection) { 
		return orthogonalMatrix;
	}
	else {
		cout<<"Orthogonal Matrix does not exist for this type of training method"<<endl;
		Array<double> empty(nBasis,nBasis);
		initArrayToZero(empty);
		return empty;
	}
}

Array<double> RBF::getOrthogonalWeightMatrix() {
	if(learner==OLS_ForwardSelection) { 
		return orthogonalWeightMatrix;
	}
	else {
		cout<<"Orthogonal Weight Matrix does not exist for this type of training method"<<endl;
		Array<double> empty(nBasis,1);
		initArrayToZero(empty);
		return empty;
	}
}

static double RBF::getDistance(Array<double> inputData, unsigned index1, Array<double> baseCentres, unsigned index2) {
	int i = inputData.dim(1);
	int j = baseCentres.dim(1);
	if(i!=j) { cout<<"error in get distance, not same dimensionality!"<<endl; exit(1);}

	int k;
	double sum=0;
	for(k=0; k<i; k++) {
		double diff=inputData(index1, k) - baseCentres(index2, k);
		sum+=diff*diff;
	}

	return sqrt(sum);
}


//======================================================================

int RBF::evaluate( Population& offsprings ) {

  unsigned InputDimension = dynamic_cast < ChromosomeT < double >& > (offsprings[0][0]).size();
  unsigned numberOfOffsprings = offsprings.size();

  Array < double > NNInput (numberOfOffsprings,InputDimension);
  Array < double > NNOutput (numberOfOffsprings, 1);
  
  for ( unsigned i = 0; i < numberOfOffsprings; i++ ) {
    for ( unsigned j = 0; j < InputDimension; j++ ) {
      NNInput(i,j) = dynamic_cast < ChromosomeT < double >& > (offsprings[i][0] )[j];
    }
  }
  evaluate(NNInput, NNOutput);

  for ( unsigned i = 0; i < numberOfOffsprings; i++ ) {
    offsprings[i].setFitness(NNOutput(i,0));
  }
  
  return 0;
}

//=======================================================================

int RBF::evaluate( Individual& offspring ) {

  unsigned InputDimension = dynamic_cast < ChromosomeT < double >& > (offspring[0]).size();
  Array < double > NNInput (1,InputDimension);
  Array < double > NNOutput (1,1);

  for ( unsigned j = 0; j < InputDimension; j++ ) {
    NNInput(0,j) = dynamic_cast < ChromosomeT < double >& > (offspring[0])[j];
  }
  evaluate(NNInput, NNOutput);

  offspring.setFitness(NNOutput(0,0));
  
  return 0;
}

int RBF::evaluate(  Array < double > inputData, Array < double >& outputData ) {
	
	unsigned nToEvaluate = inputData.dim(0);

	Array<double> localDesignMatrix(nToEvaluate, nBasis);

	unsigned i, j;


	for(i=0; i<nToEvaluate; i++) {
		for(j=0; j<nBasis; j++) {
			localDesignMatrix(i,j) = kernelFuncPtr(inputData, i, baseCentres, j, kernelParams);
		}
	}

	outputData = innerProduct(localDesignMatrix, weightMatrix);


	return 0;
}

static void RBF::KMean(unsigned nCluster, unsigned maxIter, Array<double> input, Array<double>& centresFound, Array<unsigned>& clustMap, ofstream& outFile) {
	cout<<endl<<"Start K-Mean clustering to find "<<nCluster<<" centres..."<<endl;
	outFile<<endl<<"Start K-Mean clustering to find "<<nCluster<<" centres..."<<endl;

	unsigned i, j, k;

	Rng::seed(time(NULL));
	unsigned nInput = input.dim(0);
	unsigned probDimension = input.dim(1);
	Array<double> clusterCentres(nCluster, probDimension);
	
	//randomly select one input for each cluster centre
	for(i=0; i<nCluster; i++) {
		Rng::uni.low(i*nInput/nCluster);
		Rng::uni.high((i+1)*nInput/nCluster);
		
		unsigned randNumb = (int)( floor(Rng::uni()) );
		for(j=0; j<probDimension; j++) {
			clusterCentres(i,j) = input(randNumb, j); 
		}
	}

	bool converged = false;

	Array<unsigned> clusterMap(nInput,1);
	Array<unsigned> clusterMapCopy(nInput,1);
	initArrayToZero(clusterMap);
	initArrayToZero(clusterMapCopy);

	for(i=0; i<maxIter && !converged ; i++) {
		clusterMapCopy = clusterMap;
		
		
		for(j=0; j<nInput; j++) {
			double minimDistance = getDistance(input, j, clusterCentres, clusterMap(j,0));
			for(k=0; k<nCluster ; k++) {
				if( k!=clusterMap(j,0) ) {
					double distance = getDistance(input, j, clusterCentres, k);
					if(distance<minimDistance) {
						minimDistance = distance;
						clusterMap[j] = k;
					}
				}
			}
		}

		

		//compare map
		converged = true;
		for(j=0; j<nInput; j++) {
			
			if(clusterMap(j,0)!=clusterMapCopy(j,0)) {
				converged = false;
				break;
			}
			
		}
		
		if(!converged) {
			//find new centres
			initArrayToZero(clusterCentres);
			Array<unsigned> sumArray(nCluster,1);
			initArrayToZero(sumArray);
			for(j=0; j<nInput; j++) {
				unsigned content = clusterMap(j,0);
				sumArray(content,0)++;
				for(k=0; k<probDimension; k++) {
					clusterCentres(content,k)+=input(j,k);
				}
			}
			for(j=0; j<nCluster; j++) {
				if(sumArray(j,0)>0) {
					for(k=0; k<probDimension; k++) {
						clusterCentres(j,k)/=sumArray(j,0);
					}
				}
			}
		}

		cout<<"Running K-Mean iteration "<<i<<"..."<<endl;
		outFile<<"Running K-Mean iteration "<<i<<"..."<<endl;
			
	}

	centresFound = clusterCentres;
	clustMap = clusterMap;

	cout<<"K-Mean algorithm found the following centres : "<<endl;
	cout<<"-----------------------------------------------"<<endl;

	outFile<<"K-Mean algorithm found the following centres : "<<endl;
	outFile<<"-----------------------------------------------"<<endl;

	for(i=0; i<clusterCentres.dim(0); i++) {
		cout<<"Centre "<<i<<" ==> ";
		outFile<<"Centre "<<i<<" ==> ";
		for(j=0; j<clusterCentres.dim(1); j++) {
			cout<<clusterCentres(i,j)<<"  ";
			outFile<<clusterCentres(i,j)<<"  ";
		}
		cout<<endl;
		outFile<<endl;
	}
}


