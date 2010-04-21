/*!
  \file EvoNet.h
  
  *  \author EL-Tec Group
  *  \date   Sep,2004
  *
  *
  *  \par Copyright (c) 2004: 
  *      Honda Research Institute Europe GmbH    \n
  *      Carl-Legien-Strasse 30                  \n
  *      D-63073 Offenbach/Main, Germany         \n
  *      Phone: +49 (0)69-8 90 11-758            \n
  *      Fax:   +49 (0)69-8 90 11-749            \n
  *      eMail: Stefan.Menzel@honda-ri.de        \n
  *
  *  \par Project:
  *      ApproMo
  *
  *  \par File and Revision:
  *      $RCSfile: EvoNet.h,v $
  *      $Revision: 1.9 $
  *      $Date: 2004/11/26 11:44:33 $
  *
  *
  *  \par Language and Compiler:
  *      C++, gcc (Linux)
  *
  *  \par Changes:
  *      $Log: EvoNet.h,v $
  *      Revision 1.9  2004/11/26 11:44:33  smenzel
  *      descriptions (structure optimisation) for doxygen added.
  *
  *      Revision 1.8  2004/11/26 10:09:44  smenzel
  *      added the possibility to change the parameters in the structure optimisation.
  *
  *      Revision 1.7  2004/09/15 08:10:01  smenzel
  *      the header file which fits to EvoNet.cc.
  *
  *
  *
  */
#include "MyNet.h"
#include "Population.h"
#include <string>
#include "Database.h"
#include "FFOps_OL_neu.h"

/*! \file EvoNet.h */ 
//! This class provides different functions for meta model support in EALib programs.

/*! This class includes a meta model based on Shark's Reclam and manages the neccessary database handling. For using this class in an appropriate way it is required to create an instance of the class database. After the known data is stored in the database the functionality of the neural network can be applied.*/

class EvoNet : public MyNet {
  
 public:
  
  //! dimension of input layer of the neural network
  unsigned InputDimension;
  
  //! dimension of input layer of the neural network
  unsigned OutputDimension;
  
  //! number of learning iterations
  unsigned LearnSteps;  
  
 public:
  
  /*!
    This constructor generates a network with the structure and weights given in the file named according to the string given in <Filename>.
    
    \param Filename Name of the file containing the network definition
  */
  EvoNet( string Filename );
  
  
  
  /*!
    This constructor generates a network with the parameter dimensions <in> of the input data and <out> output data and the connection matrix <cmat> and weight matrix <wmat> respectively.
    
    \param in parameter dimension of the input data
    \param out parameter dimension of the output data
    \param cmat connection matrix
    \param wmat weight matrix
  */
  EvoNet( const unsigned in, 
	  const unsigned out,
	  const Array < int > &cmat, 
	  const Array < double > &wmat );
  
  
  //! destructor
  ~EvoNet();
  
  
  //! This function optimizes the structure of the neural network based on the data stored in the arrays <Input> and <Target>. For a standard structure optimisation you can keep the default values.
  /*! 
    \param Input array containing the input data
    \param Target array containing the output data
    \param mu number of parents 
    \param lambda number of offsprings
    \param seed seed for random number generator
    \param maximum number of hidden neurons
    \param low lowest learning step size
    \param high highest learning step size
    \param sigmadel probability to remove a connection
    \param jogRate probability for a mutation of the weights
    \param B maximum number of connections connecting input and hidden neurons
    \param q tournament size
    \param number of generations
    \retval 0 ok
  */
  int StructureOptimization( Array < double > Input,
			     Array < double > Target,
			     unsigned mu = 30,
			     unsigned lambda = 30,
			     unsigned seed = 7,
			     unsigned Max_hidden = 10,
			     double low = -0.2,
			     double high = 0.2,
			     double sigmadel = 1,
			     double jogRate = 0.05,
			     unsigned B = 100,
			     unsigned q = 4,
			     unsigned Generations = 20 );
  
  //! This function optimizes the structure of the neural network based on the data stored in the database <Data>. For a standard structure optimisation you can keep the default values.
  /*! 
    \param Data the database containing the input and target data
    \param mu number of parents 
    \param lambda number of offsprings
    \param seed seed for random number generator
    \param maximum number of hidden neurons
    \param low lowest learning step size
    \param high highest learning step size
    \param sigmadel probability to remove a connection
    \param jogRate probability for a mutation of the weights
    \param B maximum number of connections connecting input and hidden neurons
    \param q tournament size
    \param number of generations
    \retval 0 ok
    \retval 1 if database <Data> is empty
  */
  int StructureOptimization( Database Data,
			     unsigned mu = 30,
			     unsigned lambda = 30,
			     unsigned seed = 7,
			     unsigned Max_hidden = 10,
			     double low = -0.2,
			     double high = 0.2,
			     double sigmadel = 1,
			     double jogRate = 0.05,
			     unsigned B = 100,
			     unsigned q = 4,
			     unsigned Generations = 20 );
  
  //! This function calculates the mean square error of the neural network considering the given database <Data>.
  /*! 
    \param Data the database containing the input and target data
    \retval MSE mean square error of the neural network
  */
  double MSE( Database Data );
  
  //! This function calculates the mean square error of the neural network considering the given data in the arrays <Input> and <Target>.
  /*! 
    \param Input array containing the input data
    \param Target array containing the target data
    \retval MSE mean square error of the neural network
  */
  double MSE( Array < double > Input,
	      Array < double > Target );
  
  //! This function calculates the mean square error of the neural network considering the given data in the Population <offsprings>.
  /*! 
    \param offsprings Population with Input data in the first Chromosome and target as the fitness value.
    \retval MSE mean square error of the neural network
  */
  double MSE( Population offsprings );
  
  
  //! This function calculates the mean square error of the neural network considering the given data in the Individual.
  /*! 
    \param offsprings Population with Input data in the first Chromosome and target as the fitness value.
    \retval MSE mean square error of the neural network
  */
  double MSE( Individual individual );
  
  //! This function evaluates the <offsprings> with the neural network and stores the results as the fitness value in each offspring.
  /*! 
    \param offsprings population of offsprings to be evaluated
    \retval 0 ok
  */
  int Evaluate( Population &offsprings );
  
  //! This function evaluates the <offspring> with the neural network and stores the results as the fitness value in the offspring.
  /*! 
    \param offspring Individual to be evaluated
    \retval 0 ok
  */
  int Evaluate( Individual &offspring );
  
  //! This function evaluates the data in <InputData> with the neural network and stores the results in <OutputData>.
  /*! 
    \param InputData Array containing the input parameters 
    \param OutputData Array containing the evaluation results
    \retval 0 ok
  */
  int Evaluate( Array < double > InputData, 
		Array < double > &OutputData );
  
  
  //! This function adapts the weights of the given neural network on the data which are stored in the database <Data>. The ReClam function Rprop is called after selecting the newest data from the database.
  
  /*!
    \param Data the data base containing the training data
    \param np    The increase factor \f$\eta^+\f$, by default set to "1.2".
    \param nm    The decrease factor \f$\eta^-\f$, by default set to "0.5".
    \param dMax  The upper limit of the increments \f$\Delta w_i^{(t)}\f$.
    The default value is "50".
    \param dMin  The lower limit of the increments \f$\Delta w_i^{(t)}\f$.
    The default value is "1e-6".
    
    \retval 0 ok
    \retval 1 no training data available
    
  */
  int Train( Database Data,
	     double np   = 1.2  ,
	     double nm   = .7   ,
	     double dMax = 50   ,
	     double dMin = 1e-6, 
	     double _delta0 = 0.01 ) 
    {
      Train_Rprop(Data,np,nm,dMax,dMin,_delta0);
      return 0;
    }
  
  
  
  //! This function adapts the weights of the given neural network on the data which are stored in the database <Data>. The ReClam function Rprop is called after selecting the newest data from the database.
  
  /*!
    \param Data the data base containing the training data
    \param np    The increase factor \f$\eta^+\f$, by default set to "1.2".
    \param nm    The decrease factor \f$\eta^-\f$, by default set to "0.5".
    \param dMax  The upper limit of the increments \f$\Delta w_i^{(t)}\f$.
    The default value is "50".
    \param dMin  The lower limit of the increments \f$\Delta w_i^{(t)}\f$.
    The default value is "1e-6".
    
    \retval 0 ok
    \retval 1 no training data available
    
  */
  
  
  
  int Train_Rprop( Database Data,
		   double np   = 1.2  ,
		   double nm   = .7   ,
		   double dMax = 50   ,
		   double dMin = 1e-6, 
		   double _delta0 = 0.01 );
  
    
    
    
    /*!
      This function adapts the weights of the given neural network on the data which are stored in the 
      \param in    The input patterns used for the training of the network.
      \param out   The target values for to the input patterns.
      \param np    The increase factor \f$\eta^+\f$, by default set to "1.2".
      \param nm    The decrease factor \f$\eta^-\f$, by default set to "0.5".
      \param dMax  The upper limit of the increments \f$\Delta w_i^{(t)}\f$. The default value is "50".
      \param dMin  The lower limit of the increments \f$\Delta w_i^{(t)}\f$. The default value is "1e-6".
      \param _delta0 Initial value for the parameter \f$\Delta\f$ The default value is "0.01".
    */
    int Train_Rprop( Array <double> InputN,
		     Array <double> TargetN,
		     int NumTrainingData,
		     double np   = 1.2  ,
		     double nm   = .7   ,
		     double dMax = 50   ,
		     double dMin = 1e-6, 
		     double _delta0 = 0.01 );

 
  //! This function scans an area of the data which must have an input dimension of 2 and an output dimension of 1. The scanning area is between the borders given as the parameters <lowerBorder> and <upperBorder>. The resolution is fixed to 100. The results can be further used for matlab.
  /*! 
    \param lowerBorder lower border of the area which is to be scanned
    \param upperBorder upper border of the area which is to be scanned
    \param OutputData Array containing the datas of the scanned areas in format [x-coordinate, y-coordinate, approximation value]
    \retval 0 ok
  */
  int ScanSquare3D ( int lowerBorder, 
		     int upperBorder, 
		     Array < double >& OutputData );

  //! This function scans an area of the data in the parameters <ind1> and <ind2> . The scanning area is between the borders given as the parameters <lowerBorder> and <upperBorder>. The resolution is fixed to 100. The results can be further used for matlab.
  /*! 
    \param lowerBorder lower border of the area which is to be scanned
    \param upperBorder upper border of the area which is to be scanned
    \param OutputData Array containing the datas of the scanned areas in format [x-coordinate, y-coordinate, approximation value]
    \param ind1 parameter 1
    \param ind2 parameter 2
    \param ParametersToKeepConstant one dimensional array containing the remaining values of the parameters which should be kept constant in increasing order
    \retval 0 ok
  */
  int ScanSquare3D( int lowerBorder,
		    int upperBorder,
		    Array < double >& OutputData,
		    int ind1,
		    int ind2,
		    Array < double > ParametersToKeepConstant );

  //! This function saves the connection and weight matrices to a file named <Filename>.
  /*! 
    \param Filename name of file in which the connection and weight matrices are stored
    \retval 0 ok
    \retval 1 if file cannot be created
  */
  int Save( string Filename );

  //! This function should loads the connection and weight matrices and generates a neural network.
  /*! 
    \param Filename name of file of which the connection and weight matrices are loaded
    \retval 0 ok
    \retval 1 if file cannot be opened
    \todo the implementation of the function
  */
  int Load( string Filename );

  //! This function returns the number of LearnSteps of the neural network.
  /*! 
    \retval LearnSteps number of LearnSteps
  */
  int GetLearnSteps();

  //! This function sets the number of LearnSteps to <NewStep>.
  /*! 
    \param NewStep new number of LearnSteps
    \retval 0 ok
  */
  int SetLearnSteps( unsigned NewStep );

  //! This function returns the dimension of the input parameters of the neural network.
  /*! 
    \retval InputDimension dimension of input parameters
  */
  int GetInputDimension();

  //! This function returns the dimension of the output parameters of the neural network.
  /*! 
    \retval OutputDimension dimension of output parameters
  */
  int GetOutputDimension();

 private:

  //! Yaochu's special Network Error (necessary in the Structure optimisation part)
  double NetworkError( Individual& Ind, 
		       FFOps_OL_neu Code, 
		       Array < double > InputData,
		       Array < double > OutputData );
 



};
