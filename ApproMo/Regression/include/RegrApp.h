#include "Database.h"
#include "Population.h"
#include "Array/Array.h"
#include "nrutil.h"
#include "ludcmp.h"
#include "lubksb.h"

 /*! \file */ 
 //! This class provides different functions for the approximation using a Regression model.
 
 /*! This class is based on the Regression approximation models. For using this class in an appropriate way it is required to create an instance of the class <database>. After the known data is stored in the database the functionality of the approximation can be applied. Therefore in the next step the model should be trained and after this the evaluation functions can be used. The consideration of the interaction terms and the degree of the Regression model can be chosen freely (>0) */

class RegrApp {

 private:

  int degOfModel;
  bool interactionTerms;
  bool InputDataIsNormalized;
  //double *ModelCoefficients;
  double *maxInColumn;
  double *minInColumn;
  //int columnDimension;
  bool ModelIsTrained;

 public:
  double *ModelCoefficients;
  int columnDimension;
  double *NaturalModelCoefficients;

  /*! \brief constructor */
  RegrApp();

  //! This constructor generates a linear regression model of order <degree> without considering interaction terms.
  /*!
    \param degree degree of the regression model (>0)
  */
  RegrApp( int degree );

  //! This constructor generates a linear regression model of order <degree>. If <interact> is TRUE the interaction terms are considered, if <FALSE> they are neglected.
  /*!
    \param degree degree of the regression model (>0)
    \param interact choice if the interaction terms shall be considered
  */
  RegrApp( int degree, 
	   bool interact );

  //! destructor
  ~RegrApp();

  //! This function calculates the mean square error of the approximation considering the given database <Data>.
  /*! 
    \param Data the database containing the input and target data
    \retval MSE mean square error of the approximation
  */
  double MSE( Database Data );

  //! This function calculates the mean square error of the approximation considering the given arrays <Input> and <Target>.
  /*! 
    \param Input array containing the input data
    \param Target array containing the target data
    \retval MSE mean square error of the approximation
  */
  double MSE( Array < double > Input,
	      Array < double > Target );

  //! This function calculates the mean square error of the approximation considering a given population <offsprings>.
  /*! 
    \param offsprings population containing the parameters in the first Chromosome and the fitness value as target value
    \retval MSE mean square error of the approximation
  */
  double MSE( Population offsprings );

  //! This function evaluates the <offsprings> with the Regression model and stores the results as the fitness value in each offspring.
  /*! 
    \param offsprings population of offsprings to be evaluated
    \retval 0 ok
    \retval 1 Model was not trained so far
  */
  int Evaluate ( Population& offsprings );

  //! This function evaluates the <offspring> with the Regression model and stores the results as the fitness value in the offspring.
  /*! 
    \param offspring Individual to be evaluated
    \retval 0 ok
    \retval 1 Model was not trained so far
  */
  int Evaluate ( Individual& offspring );

  //! This function evaluates the data in <InputData> with the Regression model and stores the results in <OutputData>.
  /*! 
    \param InputData Array containing the input parameters 
    \param OutputData Array containing the evaluation results
    \retval 0 ok
    \retval 1 Model was not trained so far
  */
  int Evaluate ( Array < double > InputData,
		 Array < double >& OutputData );

 //! This function approximates the datas which are stored in the arrays <InputData> and <TargetData> using the Regression model (default algorithm).
  /*! 
    \param InputData array containing the inputdata which are used for approximation
    \param TargetData array containing the targetdata which are used for approximation
    \retval 0 ok
    \retval 1 degree of regression model is lower 1
  */
  int Train ( Array < double > InputData,
	      Array < double > TargetData );

 //! This function approximates the datas which are stored in the arrays <InputData> and <TargetData> using the Regression model.
  /*! 
    \param InputData array containing the inputdata which are used for approximation
    \param TargetData array containing the targetdata which are used for approximation
    \retval 0 ok
    \retval 1 degree of regression model is lower 1
  */
  int Train_LSM ( Array < double > InputData,
		  Array < double > TargetData );

  //! This function approximates the datas which are stored in the database <Data> using the Regression model (default algorithm).
  /*! 
    \param Data Database containing the data which are used for approximation
    \retval 0 ok
    \retval 1 degree of regression model is lower 1
  */
  int Train ( Database Data );

  //! This function approximates the datas which are stored in the database <Data> using the Regression model.
  /*! 
    \param Data Database containing the data which are used for approximation
    \retval 0 ok
    \retval 1 degree of regression model is lower 1
  */
  int Train_LSM( Database Data );

  //! This function calculates the natural model coefficients from the trained ones (only for first order model).
  /*! 
    \param Data Database containing the data which are used for approximation
    \retval 0 ok
  */
  int CalcNaturalCoefficients();

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

  //! This function returns the degree of the regression model.
  /*! 
    \retval degOfModel degree of regression model
  */
  int GetDegreeOfRegrModel();

  //! This function sets the degree of the regression model.
  /*! 
    \param newDegree degree of the regression model
    \retval 0 
  */
  int SetDegreeOfRegrModel( int newDegree );

  //! This function returns <true> if the interaction terms are considered, <false> if not.
  /*! 
    \retval interactionTerms <true> if interaction terms are considered else <false>
  */
  bool interactionTermsConsidered();

  //! This function sets the consideration of the interaction terms of the regression model.
  /*! 
    \param choice sets the new state of the consideration of the interaction terms
    \retval 0 ok
  */
  int SetInteractionTerms( bool choice );

 private:

  int NormalizeData( Array < double > InputData,
		     Array < double >& NormInputData );

};
