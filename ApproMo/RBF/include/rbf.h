#ifndef __RBF__
#define __RBF__

#include "EALib/Population.h"
#include "Array/ArrayOp.h"




typedef enum kernel{ gaussian = 0, linear_spline, cubic_spline,  multiquadric, inverse_multiquadric, cauchy } kernelType;
typedef enum lm{ interpolate = 0, ridgeRegression, OLS_ForwardSelection } learnMethod; 
typedef enum err{ GCV = 0, UEV, FPE, BIC } errMethod;


 /*! \file */
 //! This class provides different functions for the approximation using an RBF model.

 /*! This class is based on the RBF approximation models. For using this class in an appropriate way it is required to create an instance of the class <database> or arrays of inputData and targetData. After the known data is stored in the database or the arrays, the functionality of the approximation can be applied. Therefore in the next step the model should be trained and after this the evaluation functions can be used. */


class RBF {
	private:
		unsigned inputDim;
		unsigned nBasis;
		unsigned maxNBasis;
		kernelType kernelFunction;
		errMethod errorMethod;
		learnMethod learner;
		double regularParam;	
		Array<double> kernelParams;
		double maxError;
		Array<double> weightMatrix;
		Array<double> designMatrix;
		Array<double> orthogonalMatrix;
		Array<double> orthogonalWeightMatrix;
		Array<double> baseCentres;
		unsigned maxKMeanIteration;
		unsigned maxRegIteration;
 
		double ( *kernelFuncPtr ) ( Array< double > input, unsigned i, Array< double > centre, unsigned j, Array< double > params);
		double ( *errFuncPtr ) ( Array< double > desMatrix, Array< double > output, double regParam);
		double ( *estRegParamPtr ) ( double regParam, Array< double > weightMatrix, Array< double > designMatrix, Array< double > predicted);
		void setDesignVector ( Array< double > input, unsigned index, Array< double > centres, Array< double >& outputVect);
		void setParams ( unsigned inDim, kernelType aKernel, Array< double > kernelParams);	
		
	public:
		 /*! \brief constructor */
		RBF( );		

		//! This constructor generates an Interpolation RBF model.
  		/*!
    			\param inDim dimensionality of the problem.
			\param aKernel kernel function used.
			\param kernelParams parameters for the kernel function.
  		*/
		RBF( unsigned inDim, kernelType aKernel, Array<double> kernelParams );   //interpolation

		//! This constructor generates a Ridge Regression RBF model.
		/*!
			\param inDim dimensionality of the problem.
			\param aKernel kernel function used.
			\param kernelParams parameters for the kernel function.
			\param nKernel number of kernel function centres used. 
			\param regParam starting value of the regularization parameter.
			\param maxErr maximum training MSE allowed.
			\param anErrMethod error calculation method used in optimizing the regularization parameter.
			\param maxKMeanIter maximum K-Mean clustering iteration to determine the centres of kernel functions.
			\param maxRegIter maximum regularization parameter optimization iteration.
		*/
		RBF( unsigned inDim, kernelType aKernel, Array<double> kernelParams, unsigned nKernel, double regParam, double maxErr, errMethod anErrMethod, unsigned maxKMeanIter, unsigned maxRegIter );	// ridge regression

		//! This constructor generates an Orthogonal Least Square Forward Selection RBF model.
		/*!
			\param inDim dimensionality of the problem.
			\param aKernel kernel function used.
			\param kernelParams parameters for the kernel function.
			\param maximumBasis maximum number of kernel functions.
			\param maxErr maximum training MSE allowed.
		*/
		RBF( unsigned inDim, kernelType aKernel, Array<double> kernelParams, unsigned maximumBasis,  double maxErr );		// OLS-forward selection

		 /*! \brief destructor */
		~RBF( );

		//! This function evaluates the individual <offspring> with the RBF model and stores the result as the fitness value of the offspring.
		/*!
    			\param offspring individual of offspring to be evaluated.
  		*/
		void evaluate( Individual& offspring );
		
		//! This function evaluates the population <offsprings> with the RBF model and stores the results as the fitness value of each member of offsprings.
		/*!
    			\param offsprings population of offspring to be evaluated.
  		*/
		void evaluate( Population& offsprings );
		
		//! This function evaluates the data in <InputData> with the RBF model and stores the results in <OutputData>.
		/*!
			\param inputData array containing the input parameters.
			\param outputData reference to array containing the evaluation results.
		*/
		void evaluate( Array < double > inputData, Array < double >& outputData );
		void Evaluate( Array < double > inputData, Array < double >& outputData );



		//! This function calculates the mean square error of the approximation considering the given arrays <Input> and <Target>.
		/*!
		\param Input array containing the input data.
		\param Target array containing the target data.
		\retval MSE mean square error of the approximation.
		*/
		double mse( Array < double > Input, Array < double > Target );
		double mse2( Array < double > Input, Array < double > Target );

		//! This function calculates the mean square error of the approximation considering a given population <offsprings>.
		/*!
		\param offsprings population containing the parameters in the first Chromosome and the fitness value as target value
		\retval MSE mean square error of the approximation
		*/
		double mse( Population offsprings );
		
		//! This function approximates the datas which are stored in the arrays <InputData> and <TargetData> using the RBF model.
		/*!
		\param InputData array containing the inputdata which are used for approximation.
		\param TargetData array containing the targetdata which are used for approximation.
		*/
		void train( Array < double > InputData, Array < double > TargetData );
		

		//! This function set the design(hidden) matrix of the RBF model based on <inputData> provided.
		/*!
    			\param inputData array containing the inputdata which are used for approximation.
  		*/
		void setDesignMatrix( Array< double >& inputData );
		void setDesignMatrix2( Array< double >& inputData );

		//! This function set the weight matrix of the RBF model based on <targetData> provided.
		/*!
    			\param inputData array containing the inputdata which are used for approximation.
			\param targetData array containing the target data used for approximation.
  		*/
		void setWeightMatrix( Array< double >& targetData );
		
		void setWeightMatrix2( Array< double >& targetData );


		//! This function returns the design(hidden) matrix of the RBF model.
		/*!
    			\retval designMatrix array containing the <designMatrix> of the RBF model.
  		*/
		Array<double> getDesignMatrix( );

		//! This function returns the weight matrix of the RBF model.
		/*!
    			\retval weightMatrix array containing the <weightMatrix> of the RBF model.
  		*/
		Array<double> getWeightMatrix( );

		//! This function returns the orthogonal matrix of the OLSForward (Orthogonal Least Square - Forward Selection) RBF model.
		/*!
    			\retval orthogonalMatrix array containing the <orthogonalMatrix> of the RBF model.
  		*/
		Array<double> getOrthogonalMatrix( );

		//! This function returns the orthogonal weight matrix of the OLSForward (Orthogonal Least Square - Forward Selection) RBF model.
		/*!
    			\retval orthogonalWeightMatrix array containing the <orthogonalWeightMatrix> of the RBF model.
  		*/
		Array<double> getOrthogonalWeightMatrix( );
		
		//! This function returns the kernel function centres of the RBF model.
		/*!
    			\retval baseCentres array containing the <baseCentres> of the RBF model.
  		*/
		Array<double> getBaseCentres( );

		//! This function calculates the distance between <i>-th <input> and j-th <centre>.
		/*!
			\param input array containing the input which are used for training.
			\param i index of the input which distance to be calculated.
			\param centre array containing the centres of the RBF model.
			\param j index of the centre to which the distance from input is to be calculated.
  			\retval distance value of distance.
		*/
		static double getDistance( Array< double >& input, unsigned i, Array<double>& centre, unsigned j );

		//! This function calculates the gaussian kernel output between <i>-th <input> and j-th <centre>.
		/*!
			\param input array containing the input which are used for training.
			\param i index of the input which distance to be calculated.
			\param centre array containing the centres of the RBF model.
			\param j index of the centre to which the distance from input is to be calculated.
			\param params array containing parameter(s) for the gaussian kernel function.
  			\retval kernelOutput value of kernel output.
		*/
		static double gaussFunc( Array< double > input, unsigned i, Array< double > centre, unsigned j, Array< double > params );

		//! This function calculates the cubic spline kernel output between <i>-th <input> and j-th <centre>.
		/*!
			\param input array containing the input which are used for training.
			\param i index of the input which distance to be calculated.
			\param centre array containing the centres of the RBF model.
			\param j index of the centre to which the distance from input is to be calculated.
			\param params array containing parameter(s) for the cubic spline kernel function.
  			\retval kernelOutput value of kernel output.
		*/
		static double cubicFunc( Array< double > input, unsigned i, Array< double > centre, unsigned j, Array< double > params );

		//! This function calculates the linear spline kernel output between <i>-th <input> and j-th <centre>.
		/*!
			\param input array containing the input which are used for training.
			\param i index of the input which distance to be calculated.
			\param centre array containing the centres of the RBF model.
			\param j index of the centre to which the distance from input is to be calculated.
			\param params array containing parameter(s) for the linear spline kernel function.
  			\retval kernelOutput value of kernel output.
		*/
		static double linearFunc( Array< double > input, unsigned i, Array< double > centre, unsigned j, Array< double > params );


		//! This function calculates the multiquadrics kernel output between <i>-th <input> and j-th <centre>.
		/*!
			\param input array containing the input which are used for training.
			\param i index of the input which distance to be calculated.
			\param centre array containing the centres of the RBF model.
			\param j index of the centre to which the distance from input is to be calculated.
			\param params array containing parameter(s) for the multiquadrics kernel function.
  			\retval kernelOutput value of kernel output.
		*/
		static double multiquadricFunc( Array< double > input, unsigned i, Array< double > centre, unsigned j, Array< double > params);

		//! This function calculates the inverse multiquadrics kernel output between <i>-th <input> and j-th <centre>.
		/*!
			\param input array containing the input which are used for training.
			\param i index of the input which distance to be calculated.
			\param centre array containing the centres of the RBF model.
			\param j index of the centre to which the distance from input is to be calculated.
			\param params array containing parameter(s) for the inverse multiquadrics kernel function.
  			\retval kernelOutput value of kernel output.
		*/
		static double inverseMultiquadricFunc( Array< double > input, unsigned i, Array< double > centre, unsigned j, Array< double > params);

		//! This function calculates the cauchy kernel output between <i>-th <input> and j-th <centre>.
		/*!
			\param input array containing the input which are used for training.
			\param i index of the input which distance to be calculated.
			\param centre array containing the centres of the RBF model.
			\param j index of the centre to which the distance from input is to be calculated.
			\param params array containing parameter(s) for the cauchy kernel function.
			\retval kernelOutput value of kernel output.
  		*/
		static double cauchyFunc( Array< double > input, unsigned i, Array< double > centre, unsigned j, Array< double > params );		

		//! This function is a helper function to calculate sum squared error analytically. It is used in the regularization parameter optimization.
		/*!
			\param desMatrix array containing the design matrix.
			\param output array containing the desired output value.
			\param regParam regularization parameter.
			\param A reference to array containing a helper matrix A.	
  			\retval sumSquareError sum squared error
		*/
		static double YPYFunc( Array< double > desMatrix, Array< double > output, double regParam, Array< double >& A ); //Compute Sum Squared Error analytically

		//! This function calculates the Generalized Cross-Validation error.
		/*!
			\param desMatrix array containing the design matrix.
			\param output array containing the desired output value.
			\param regParam regularization parameter.
			\retval GCVError GCV error.
  		*/
		static double GCVFunc( Array< double > desMatrix, Array< double > output, double regParam ); // Generalized Cross-Validation

		//! This function calculates the Unbiased Estimate Variance.
		/*!
			\param desMatrix array containing the design matrix.
			\param output array containing the desired output value.
			\param regParam regularization parameter.
			\retval UEVError UEV error.
  		*/
		static double UEVFunc( Array< double > desMatrix, Array< double > output, double regParam ); // Unbiased Estimate Variance

		//! This function calculates the Final Prediction Error.
		/*!
			\param desMatrix array containing the design matrix.
			\param output array containing the desired output value.
			\param regParam regularization parameter.
			\retval FPEError FPE error.
  		*/
		static double FPEFunc( Array< double > desMatrix, Array< double > output, double regParam ); // Final Prediction Error
		
		//! This function calculates the Bayesian Information Criterion error.
		/*!
			\param desMatrix array containing the design matrix.
			\param output array containing the desired output value.
			\param regParam regularization parameter.
			\retval BICError BIC error.
  		*/
		static double BICFunc( Array< double > desMatrix, Array< double > output, double regParam ); // Bayesian Information Criterion

		//! This function estimates the new regularization parameter based on the Generalized Cross-Validation error.
		/*!
			\param regParam regularization parameter.
			\param weightMatrix array containing the weight matrix.
			\param designMatrix array containing the design matrix.
			\param predicted array containing the predicted output.
			\retval newRegParam new regularization parameter.
  		*/
		static double estRegParamGCV( double regParam, Array< double > weightMatrix, Array< double > designMatrix, Array< double > predicted);
		
		//! This function estimates the new regularization parameter based on the Unbiased Estimate Variance.
		/*!
			\param regParam regularization parameter.
			\param weightMatrix array containing the weight matrix.
			\param designMatrix array containing the design matrix.
			\param predicted array containing the predicted output.
			\retval newRegParam new regularization parameter.
  		*/
		static double estRegParamUEV( double regParam, Array< double > weightMatrix, Array< double > designMatrix, Array< double > predicted);

		//! This function estimates the new regularization parameter based on the Final Prediction Error.
		/*!
			\param regParam regularization parameter.
			\param weightMatrix array containing the weight matrix.
			\param designMatrix array containing the design matrix.
			\param predicted array containing the predicted output.
			\retval newRegParam new regularization parameter.
  		*/
		static double estRegParamFPE( double regParam, Array< double > weightMatrix, Array< double > designMatrix, Array< double > predicted);

		//! This function estimates the new regularization parameter based on the Bayesian Information Criterion.
		/*!
			\param regParam regularization parameter.
			\param weightMatrix array containing the weight matrix.
			\param designMatrix array containing the design matrix.
			\param predicted array containing the predicted output.
			\retval newRegParam new regularization parameter.
  		*/
		static double estRegParamBIC( double regParam, Array< double > weightMatrix, Array< double > designMatrix, Array< double > predicted);

		//! This function performs the K-Mean clustering algorithm to provide <nCluster> clusters based on a set of inputs <input>.
		/*!
			\param nCluster number of desired clusters.
			\param maxIter maximum iteration count.
			\param input array containing the inputs.
			\param centresFound reference to array containing the cluster centres found.
			\param clustMap reference to array containing the pair of inputs and cluster numbers to which each of them are assigned.
  		*/
		static void KMean( unsigned nCluster, unsigned maxIter, Array< double > input, Array< double >& centresFound, Array< unsigned >& clustMap);
		
		//! This function initializes an array to zero.
		/*!
			\param arr array to be initialized to zero.
  		*/
		template <class Type>
		static void initArrayToZero( Array< Type >& arr ) {
		
			unsigned i, j;
			for ( i = 0; i < arr.dim( 0 ); i++ ) {
				for ( j = 0; j < arr.dim( 1 ); j++ ) {
					arr( i, j ) = 0;
				}
			}	
		}
};

#endif
