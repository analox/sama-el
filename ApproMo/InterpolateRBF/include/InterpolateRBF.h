#ifndef __INTERPOLATE_RBF__
#define __INTERPOLATE_RBF__

#include "EALib/Population.h"
#include "Array/ArrayOp.h"


typedef enum kernel{ gaussian = 0, linear_spline, cubic_spline,  multiquadric, inverse_multiquadric, cauchy } kernelType;

 /*! \file */
 //! This class provides different functions for the approximation using an InterpolationRBF model.

 /*! This class is based on the RBF approximation models.  Therefore in the next step the model should be trained and after this the evaluation functions can be used. */


class InterpolateRBF {
	private:
		unsigned inputDim;
		unsigned nBasis;
		unsigned maxNBasis;
		kernelType kernelFunction;
		Array<double> kernelParams;
		double maxError;
		Array<double> weightMatrix;
		Array<double> designMatrix;
		Array<double> baseCentres;
 
		double ( *kernelFuncPtr ) ( Array< double > input, unsigned i, Array< double > centre, unsigned j, Array< double > params);
		void setDesignVector ( Array< double > input, unsigned index, Array< double > centres, Array< double >& outputVect);
		void setParams ( unsigned inDim, kernelType aKernel, Array< double > kernelParams);	
		
	public:
		 /*! \brief constructor */
		InterpolateRBF( );		

		//! This constructor generates an Interpolation RBF model.
  		/*!
    			\param inDim dimensionality of the problem.
			\param aKernel kernel function used.
			\param kernelParams parameters for the kernel function.
  		*/
		InterpolateRBF( unsigned inDim, kernelType aKernel, Array<double> kernelParams );   //interpolation

		 /*! \brief destructor */
		~InterpolateRBF( );

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



		//! This function calculates the mean square error of the approximation considering the given arrays <Input> and <Target>.
		/*!
		\param Input array containing the input data.
		\param Target array containing the target data.
		\retval MSE mean square error of the approximation.
		*/
		double mse( Array < double > Input, Array < double > Target );

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

		//! This function set the weight matrix of the RBF model based on <targetData> provided.
		/*!
    			\param inputData array containing the inputdata which are used for approximation.
			\param targetData array containing the target data used for approximation.
  		*/
		void setWeightMatrix( Array< double >& targetData );
		


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
		double getDistance( Array< double >& input, unsigned i, Array<double>& centre, unsigned j );

		//! This function calculates the gaussian kernel output between <i>-th <input> and j-th <centre>.
		/*!
			\param input array containing the input which are used for training.
			\param i index of the input which distance to be calculated.
			\param centre array containing the centres of the RBF model.
			\param j index of the centre to which the distance from input is to be calculated.
			\param params array containing parameter(s) for the gaussian kernel function.
  			\retval kernelOutput value of kernel output.
		*/
		double gaussFunc( Array< double > input, unsigned i, Array< double > centre, unsigned j, Array< double > params );

		//! This function calculates the cubic spline kernel output between <i>-th <input> and j-th <centre>.
		/*!
			\param input array containing the input which are used for training.
			\param i index of the input which distance to be calculated.
			\param centre array containing the centres of the RBF model.
			\param j index of the centre to which the distance from input is to be calculated.
			\param params array containing parameter(s) for the cubic spline kernel function.
  			\retval kernelOutput value of kernel output.
		*/
		double cubicFunc( Array< double > input, unsigned i, Array< double > centre, unsigned j, Array< double > params );

		//! This function calculates the linear spline kernel output between <i>-th <input> and j-th <centre>.
		/*!
			\param input array containing the input which are used for training.
			\param i index of the input which distance to be calculated.
			\param centre array containing the centres of the RBF model.
			\param j index of the centre to which the distance from input is to be calculated.
			\param params array containing parameter(s) for the linear spline kernel function.
  			\retval kernelOutput value of kernel output.
		*/
		double linearFunc( Array< double > input, unsigned i, Array< double > centre, unsigned j, Array< double > params );


		//! This function calculates the multiquadrics kernel output between <i>-th <input> and j-th <centre>.
		/*!
			\param input array containing the input which are used for training.
			\param i index of the input which distance to be calculated.
			\param centre array containing the centres of the RBF model.
			\param j index of the centre to which the distance from input is to be calculated.
			\param params array containing parameter(s) for the multiquadrics kernel function.
  			\retval kernelOutput value of kernel output.
		*/
		double multiquadricFunc( Array< double > input, unsigned i, Array< double > centre, unsigned j, Array< double > params);

		//! This function calculates the inverse multiquadrics kernel output between <i>-th <input> and j-th <centre>.
		/*!
			\param input array containing the input which are used for training.
			\param i index of the input which distance to be calculated.
			\param centre array containing the centres of the RBF model.
			\param j index of the centre to which the distance from input is to be calculated.
			\param params array containing parameter(s) for the inverse multiquadrics kernel function.
  			\retval kernelOutput value of kernel output.
		*/
		double inverseMultiquadricFunc( Array< double > input, unsigned i, Array< double > centre, unsigned j, Array< double > params);

		//! This function calculates the cauchy kernel output between <i>-th <input> and j-th <centre>.
		/*!
			\param input array containing the input which are used for training.
			\param i index of the input which distance to be calculated.
			\param centre array containing the centres of the RBF model.
			\param j index of the centre to which the distance from input is to be calculated.
			\param params array containing parameter(s) for the cauchy kernel function.
			\retval kernelOutput value of kernel output.
  		*/
		double cauchyFunc( Array< double > input, unsigned i, Array< double > centre, unsigned j, Array< double > params );		


		//! This function performs the K-Mean clustering algorithm to provide <nCluster> clusters based on a set of inputs <input>.
		/*!
			\param nCluster number of desired clusters.
			\param maxIter maximum iteration count.
			\param input array containing the inputs.
			\param centresFound reference to array containing the cluster centres found.
			\param clustMap reference to array containing the pair of inputs and cluster numbers to which each of them are assigned.
  		*/
		void KMean( unsigned nCluster, unsigned maxIter, Array< double > input, Array< double >& centresFound, Array< unsigned >& clustMap);
		
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
