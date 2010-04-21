/*!
  \file EvoNet.cc
  
  *  \author EL-Tec Group
  *  \date   Sep,2004
  *
  *
  *  \par Copyright (c) 2004: 
  *      Honda Research Institute Europe GmbH    \n
  *      Carl-Legien-Strasse 30                  \n
  *      D-63073 Offenbach/Main, Germany         \n
  *      Phone: +49 (0)69-8 90 11-734            \n
  *      Fax:   +49 (0)69-8 90 11-749            \n
  *      eMail: Stefan.Menzel@honda-ri.de       \n
  *
  *  \par Project:
  *      ApproMo
  *
  *  \par File and Revision:
  *      $RCSfile: EvoNet.cc,v $
  *      $Revision: 1.13 $
  *      $Date: 2004/11/26 10:09:44 $
  *
  *
  *  \par Language and Compiler:
  *      C++, gcc (Linux)
  *
  *  \par Changes:
  *      $Log: EvoNet.cc,v $
  *      Revision 1.13  2004/11/26 10:09:44  smenzel
  *      added the possibility to change the parameters in the structure optimisation.
  *
  *      Revision 1.12  2004/09/15 08:22:30  smenzel
  *      Warnings revised.
  *
  *      Revision 1.11  2004/09/14 16:21:14  smenzel
  *      the most recent version of EvoNet.
  *
  *      Revision 1.10  2004/09/14 14:59:19  markus
  *      header for structure optimisation changed
  *      parameters are now defined in the source code
  *
  *      Revision 1.9  2004/09/14 14:51:48  markus
  *      now after the update the version look somehow strange...
  *      maybe we are on the way going back one or two versions...
  *      anyway, a new version with two errors less
  *
  *
  *
  */
#include "EvoNet.h"
#include "sqr.h"

#define LAMARCK //Use lamarckian mechanism, otherwise, Baldwin effect 

#undef DEBUG

using namespace std;

//=======================================================================

EvoNet::EvoNet( string Filename ) : MyNet ( Filename ) {
  
  LearnSteps = 100;

  OutputDimension = outputDimension;
  InputDimension  = inputDimension;
}

//=======================================================================

EvoNet::EvoNet( const unsigned in, 
		const unsigned out,
		const Array < int > &cmat, 
		const Array < double > &wmat ) : MyNet ( in, out, cmat, wmat ) {
 
  LearnSteps = 100;
    
  OutputDimension = out;
  InputDimension = in;
}

//=======================================================================

EvoNet::~EvoNet() {
}

//=======================================================================

int EvoNet::StructureOptimization ( Array < double > Input,
				    Array < double > Target,
				    unsigned mu,
				    unsigned lambda,
				    unsigned seed,
				    unsigned Max_hidden,
				    double low,
				    double high,
				    double sigmadel,
				    double jogRate,
				    unsigned B,
				    unsigned q,
				    unsigned Generations ) {

  Database tempData;
  tempData.AddTrainingData ( Input, Target );
  StructureOptimization ( tempData, 
			  mu,
			  lambda,
			  seed,
			  Max_hidden,
			  low,
			  high,
			  sigmadel,
			  jogRate,
			  B,
			  q,
			  Generations );
  return 0;
}

//=======================================================================


int EvoNet::StructureOptimization( Database Data,
				   unsigned mu,
				   unsigned lambda,
				   unsigned seed,
				   unsigned Max_hidden,
				   double low,
				   double high,
				   double sigmadel,
				   double jogRate,
				   unsigned B,
				   unsigned q,
				   unsigned Generations ) {
  
  cout << "Structure Optimisation " << endl;
    
  // ---- prepare training data (copy them out of the data base)
  //              (only the newest training data are used)
    
  if (Data.GetNumTrainingData() == 0) {
    cout << "No training data available!" << endl;
    return 1;
  }

  int UsedTrainingData = Data.GetUsedTrainingData();
  int DataInDatabase = Data.GetNumTrainingData();
  
  cout << "Used Training Data:      " << UsedTrainingData << endl;
  cout << "Available Training Data: " << DataInDatabase << endl;
  
  Array < double > InputData(UsedTrainingData, InputDimension);
  Array < double > OutputData(UsedTrainingData, OutputDimension);
 
  Data.GetLastData(InputData,OutputData);
  
#ifdef DEBUG
    cout << "Training Data" << endl;
    cout << "Input Dimension: " << InputDimension << "Output Dimension: " << OutputDimension << endl;
    for (int i = 0; i < InputData.dim(0); i++) {
      for (unsigned int j = 0; j < InputDimension; j++)
	cout << InputData(i,j) << " ";
      cout << "|" << OutputData(i,0) << endl;
    }
#endif

  // start Yaochu's magic
  
  Individual Ind_best;
  FFOps_OL_neu Code;
  
  // define parameters for learning
  /*  unsigned mu = 30;
  unsigned lambda = 30;
  unsigned seed = 7;
  unsigned Max_hidden = 10;                // Specify the maximal number of hidden nodes
  */
  unsigned N_input = InputDimension;       // Input dimension
  unsigned N_output = OutputDimension;     // output dimension
  /*
  double low = -0.2;                       // maximal negative step-size
  double high = 0.2;                       // minimal positive step-size
  double sigmadel = 1;
  double jogRate = 0.05;
 
  unsigned B = 100;
  unsigned q = 4;
  unsigned Generations = 20;
  */
  // cout << "Initialize population." << endl;

  Population parents, offsprings;
  parents   .resize(mu);
  parents   .setMinimize( );
  offsprings.resize(lambda);
  offsprings.setMinimize( ); 

  Rng::seed(seed);

  for( unsigned i = 0; i < mu; i++ ){
    int units=Rng::discrete(1,Max_hidden);
    parents[i] = Code.createEmpty(N_input, units, N_output);
    double A=(double)(B-units*N_output-units-N_output)/(double)(N_input*units);
    Code.initConnections(parents[i],A,1,1);
    Code.initWeights(parents[i],low,high);
  }
  // evaluate their initial fitness
  for( unsigned i = 0; i < mu; i++ ){
    double FitN = NetworkError(parents[i], Code, InputData, OutputData);
    parents[i].setFitness(FitN);
  }
 
  ///////////////////////////////////////////////////////////////////////
  //
  // iterate
  //
  ///////////////////////////////////////////////////////////////////////
  
  //cout << "Iteration begins." << endl;
  for( unsigned t = 1; t <= Generations; t++ ) {

    //reproduce
    for (unsigned i=0; i<lambda;i++)
      offsprings[i] = parents[i];
    
    // mutate
    for (unsigned i=0; i < lambda; i++){
      unsigned OpIndex = Rng::discrete(0,4);
      switch(OpIndex){
      case 0:
	Code.addConnection(offsprings[i], low, high);
	break; 
      case 1:
	Code.deleteConnection(offsprings[i],sigmadel);
	break; 
      case 2:
	Code.jogWeights(offsprings[i],jogRate);
	break;
      case 3:
	if (Code.getNHid(offsprings[i]) < Max_hidden)
	  Code.addNeuron(offsprings[i],low, high);
	break;
      case 4:
	Code.deleteNeuron(offsprings[i]);
	break;
      }
    }
      
    // evaluate offspring's fitness
    for( unsigned i = 0; i < lambda; i++ ){
      //   MARKUS: Windows has a compiler error here 
      double FitN = NetworkError(offsprings[i], Code, InputData, OutputData);
      offsprings[i].setFitness(FitN);
    }

    // selection
    //
    parents.selectEPTournament(offsprings,(unsigned)q);
    
    // cout << t << "   " << parents.minFitness() << endl;
  }
  Ind_best = parents.best();
  
  N_input  = Code.getNIn(Ind_best);
  N_output = Code.getNOut(Ind_best);
  Array<int>    con = Code.getConnectionsReClaM(Ind_best);
  Array<double> wei = Code.getWeightsReClaM(Ind_best);
  
  setStructure ( con, wei ); 
  cout << "size of wei " << wei.dim(0)<< endl;
  cout << "size of con " << con.dim(0)<< endl;
  
  return 0;
}

//=======================================================================

double EvoNet::MSE( Database Data ) {
  
  int UsedTrainingData = Data.GetUsedTrainingData();
  Array < double > InputData(UsedTrainingData, InputDimension);
  Array < double > TargetData(UsedTrainingData, OutputDimension);
 
  Data.GetLastData(InputData,TargetData);
  
  return MSE(InputData,TargetData);
}

//=======================================================================


double EvoNet::MSE( Population offsprings ) {
  
  // input data are expected in the first chromosome
  unsigned NNDimension = dynamic_cast< ChromosomeT< double >& >(offsprings[0][0]).size();
  unsigned numberOfOffsprings = offsprings.size();

  Array < double > NNInput(numberOfOffsprings, InputDimension);
  Array < double > NNTarget(numberOfOffsprings, OutputDimension);
 
  for( unsigned i = 0; i < offsprings.size( ); i++ ) {  
    for(unsigned jj = 0; jj < NNDimension; jj++) {
      NNInput(i,jj) = dynamic_cast< ChromosomeT< double >& >( offsprings[i][0] )[jj];
    }
    NNTarget(i,0) = offsprings[i].fitnessValue();
  }

  return MSE(NNInput,NNTarget);
}
//=======================================================================

double EvoNet::MSE( Individual individual ) {
  
  // input data are expected in the first chromosome
  unsigned NNDimension = dynamic_cast< ChromosomeT< double >& >(individual[0]).size();
  unsigned numberOfOffsprings = 1;
  
  Array < double > NNInput(numberOfOffsprings, InputDimension);
  Array < double > NNTarget(numberOfOffsprings, OutputDimension);
  
  for(unsigned jj = 0; jj < NNDimension; jj++) {
    NNInput(0,jj) = dynamic_cast< ChromosomeT< double >& >( individual[0] )[jj];
  }
  NNTarget(0,0) = individual.fitnessValue();

  return MSE(NNInput,NNTarget);
}

//=======================================================================

double EvoNet::MSE( Array < double > InputN,
		    Array < double > TargetN ) {
  
  return meanSquaredError(InputN,TargetN); 
}

//=======================================================================

int EvoNet::Evaluate( Population &offsprings ) {

// input data are expected in the first chromosome
unsigned NNDimension = dynamic_cast< ChromosomeT< double >& >(offsprings[0][0]).size();
  
  Array < double > NNInput(InputDimension);
  Array < double > NNOutput(OutputDimension);
 
  // --- Use the network for all individuals
  for( unsigned i = 0; i < offsprings.size( ); i++ ) {  
    NNInput.resize(NNDimension); 
    
    for(unsigned jj = 0; jj < NNDimension; jj++) 
      NNInput[jj] = dynamic_cast< ChromosomeT< double >& >( offsprings[i][0] )[jj];

    model(NNInput, NNOutput);
    double fitness = NNOutput(0);

    offsprings[i].setFitness(fitness);
  }

  return 0;
}

//=======================================================================

int EvoNet::Evaluate( Individual &offspring ) {
  
  // training data will be the first chromosome
  int NNDimension = dynamic_cast< ChromosomeT< double >& >(offspring[0]).size();
  
  Array < double > NNInput(InputDimension);
  Array < double > NNOutput(OutputDimension);
 
  NNInput.resize(NNDimension); 
    
  for(int jj = 0; jj < NNDimension; jj++) 
    NNInput[jj] = dynamic_cast< ChromosomeT< double >& >( offspring[0] )[jj];
  
  model(NNInput, NNOutput);
  double fitness = NNOutput(0);

  offspring.setFitness(fitness);

  return 0;
}

//=======================================================================

int EvoNet::Evaluate( Array < double > InputData, 
		      Array < double > &OutputData ) {
  
  Array < double > NNInput ( InputDimension );
  Array < double > NNOutput ( OutputDimension );
  OutputData.resize ( (unsigned int)InputData.dim(0), OutputDimension );
  for ( int i = 0; i < (int)InputData.dim(0); i++ ) {
    for(int jj = 0; jj < (int)InputDimension; jj++) {
      NNInput[jj] = InputData ( i,jj );
    }
    model(NNInput, NNOutput);
    for(int jj = 0; jj < (int)OutputDimension; jj++) {
      OutputData(i,jj) = NNOutput(jj);
    }
  }
  //  model(InputData, OutputData);
  
  return 0;
}

//=======================================================================

int EvoNet::ScanSquare3D( int lowerBorder,
			  int upperBorder,
			  Array < double >& OutputData ) {

  int i = 0;
  int j = 1;
  Array < double > values (0);
  ScanSquare3D ( lowerBorder,
		 upperBorder,
		 OutputData,
		 i,
		 j,
		 values );
  return 0;
}

//=======================================================================

int EvoNet::ScanSquare3D( int lowerBorder,
			  int upperBorder,
			  Array < double >& OutputData,
			  int ind1,
			  int ind2,
			  Array < double > ParametersToKeepConstant ) {

  Array < double > InputData( 1,2 );
  Array < double > Output( 1,1 );
  Array < double > DataToEvaluate ( 1, ParametersToKeepConstant.dim(0)+2 );

  int resolution = 100;
  int diff = upperBorder - lowerBorder;
  OutputData.resize( (unsigned int)sqr(resolution+1), (unsigned int)3);
  int outputEntry = 0;
  for ( int i = 0; i <= resolution; i++ ) {
    for ( int j = 0; j <= resolution; j++ ) {
      InputData(0,0) = (double)i*diff/resolution+lowerBorder;
      InputData(0,1) = (double)j*diff/resolution+lowerBorder;
      OutputData(outputEntry,0) = (double)i*diff/resolution+lowerBorder;
      OutputData(outputEntry,1) = (double)j*diff/resolution+lowerBorder;
      
      if ( DataToEvaluate.dim(1) > 2 ) {
	int counter = 0;
	for ( int k = 0; k < (int)DataToEvaluate.dim(1); k++ ) {
	  if ( k == ind1 ) {
	    DataToEvaluate ( 0, k ) = InputData (0,0);
	  }
	  else {
	    if ( k == ind2 ) {
	      DataToEvaluate ( 0, k ) = InputData ( 0, 1 );
	    }
	    else {
	      DataToEvaluate ( 0, k ) = ParametersToKeepConstant ( counter );
	      counter++;
	    }
	  }
	}
      }
      else {
	DataToEvaluate ( 0, 0 ) = InputData ( 0, 0 );
	DataToEvaluate ( 0, 1 ) = InputData ( 0, 1 );
      }

      Evaluate(DataToEvaluate,Output);
      OutputData(outputEntry,2) = Output(0,0);
      outputEntry++;
    }
  }

  return 0;
}

//=======================================================================

int EvoNet::Train_Rprop( Database Data,
			 double np,
			 double nm,
			 double dMax,
			 double dMin, 
			 double _delta0 ) 
{
  
  if (Data.GetUsedTrainingData() == 0) {
    cout << "No training data available!" << endl;
    return 1;
  }
  
  int usedData  = Data.GetUsedTrainingData();
  int totalData = Data.GetNumTrainingData();
  
  Array < double > InputData(usedData, InputDimension);
  Array < double > TargetData(usedData, OutputDimension);
  
  Data.GetLastData(InputData, TargetData);
  
  Train_Rprop(InputData, TargetData, totalData, np, nm, dMax, dMin, _delta0 );
  
  return 0;
}

//=======================================================================

int EvoNet::Train_Rprop( Array <double> InputN, 
			 Array <double> TargetN,
			 int NumTrainingData,
			 double np,
			 double nm,
			 double dMax,
			 double dMin, 
			 double _delta0) 
{
  

#ifdef DEBUG
  cout << "Number of Training Data: " << NumTrainingData << endl;
  cout << "Number of Used Training Data: " << InputN.dim(0) << endl;
  cout << "Number of Learn Steps:   " << LearnSteps << endl;
#endif
  
#ifdef DEBUG  
  Array <double> OutputN (InputN.dim(0), OutputDimension);
  model(InputN, OutputN);
  double rms_LB = 0.;
  for (unsigned i = 0; i < InputN.dim(0); i++){
    rms_LB += (OutputN(i,0) - TargetN(i,0)) * (OutputN(i,0) - TargetN(i,0));
  }
  printf("Before training: %f ", sqrt(rms_LB/(double)NumTrainingData));
#endif
  
  initRprop(_delta0);
  
  for (unsigned iteration = 0; iteration < (unsigned)LearnSteps; iteration++)
    rprop(InputN, TargetN, np, nm, dMax, dMin);
  
#ifdef DEBUG
  model(InputN,OutputN);
  double rms_LA = 0.;
  for (unsigned i = 0; i < InputN.dim(0); i++){
    rms_LA += (OutputN(i,0) - TargetN(i,0)) * (OutputN(i,0) - TargetN(i,0));
  }
  
  printf("After Training:%f \n", sqrt(rms_LA/(double)NumTrainingData ));
#endif  
  return 0;
}

//=======================================================================

int EvoNet::Save( string Filename ) {

  // write the network
  FILE* fp = fopen(Filename.c_str(),"w");
  if (fp == NULL) return 1;
  fprintf(fp,"%d %d\n",inputDimension,outputDimension);
  fprintf(fp,"\n");
  
  Array <int> c=getConnections();
  
  for (unsigned i=0; i<c.dim(0); i++){
    for (unsigned j=0; j<c.dim(1)-1; j++)
      fprintf(fp,"%d ",c(i,j));
    fprintf(fp,"%d\n",c(i,c.dim(1)-1));
  }
  
  fprintf(fp,"\n");
  Array<double> w=getWeights();
  
  for (unsigned i=0; i<w.dim(0); i++){
    for (unsigned j=0; j<w.dim(1)-1; j++)
      fprintf(fp,"%10.8e ",w(i,j));
    fprintf(fp,"%10.8e\n",w(i,w.dim(1)-1));
  }
  fflush(fp);
  fclose(fp);
  
  return 0;
}

//=======================================================================

int EvoNet::Load( string Filename ) {

  /*
  metaNet = new MyNet(Filename);
  LearnSteps = 100;
  
  OutputDimension = metaNet->outputDimension;
  InputDimension  = metaNet->inputDimension;
  */
  return 0;
}

//=======================================================================

int EvoNet::GetLearnSteps() {

  return LearnSteps;

}

//=======================================================================

int EvoNet::SetLearnSteps( unsigned NewStep ) {

  LearnSteps = NewStep;

  return 0;
}

//=======================================================================

int EvoNet::GetInputDimension() {

  return InputDimension;

}

//=======================================================================

int EvoNet::GetOutputDimension() {

  return OutputDimension;

}



//======================================================================= 
//     this function is needed in the Structure Optimisation
//======================================================================= 

double EvoNet::NetworkError( Individual& Ind, 
				   FFOps_OL_neu Code, 
				   Array < double > InputData,
				   Array < double > OutputData ) {

  double np = 1.2;
  double nm = 0.5;
  int cycles = 50;
  double delta0 = 0.01;
  double dMin = 0;       // minimal learning rate in rprop
  double dMax = 50.;     // maximal learning rate in the rprop 
  double error;

  // initialize variables
  unsigned N_input  = Code.getNIn(Ind);
  unsigned N_output = Code.getNOut(Ind);
  Array<int>    con = Code.getConnectionsReClaM(Ind);
  Array<double> wei = Code.getWeightsReClaM(Ind);

  MyNet net(N_input,N_output,con,wei);
  
  net.initRprop(delta0);
  for (int iteration=0; iteration<cycles; iteration++)
    net.rprop(InputData, OutputData, np, nm, dMax, dMin);
 
  error = net.meanSquaredError(InputData, OutputData);

#ifdef LAMARCK
  Code.setWeightsReClaM(Ind, net.getWeights());
#endif
  return error;
}
