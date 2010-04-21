#include "Database.h"
#include "Population.h"
#include "ArrayOp.h"
#include <fstream.h>

typedef enum kernel{ gaussian=0, linear_spline, cubic_spline, thin_plate_spline, multiquadric, inverse_multiquadric, cauchy} kernelType;
typedef enum lm{ interpolate=0, ridgeRegression, OLS_ForwardSelection} learnMethod; 
typedef enum err{ GCV=0, UEV, FPE, BIC } errMethod;

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
		double (*kernelFuncPtr)(Array<double> input, unsigned i, Array<double> centre, unsigned j, Array<double> params);
		double (*errFuncPtr)(Array<double> desMatrix, Array<double> output, double regParam);
		double (*estRegParamPtr)(double regParam, Array<double> weightMatrix, Array<double> designMatrix, Array<double>predicted);
		void setDesignVector(Array<double> input, unsigned index, Array<double> centres, Array<double>& outputVect);
		void setParams(unsigned inDim, kernelType aKernel, Array<double> kernelParams);
		
		static void printToFile(unsigned iter, Array<double> wMatrix, Array<double> centres, Array<double> kernParams, kernelType kernFunction, learnMethod aMethod);
		

	public:
		RBF( );		//use default setting;
		RBF( unsigned inDim, kernelType aKernel, Array<double> kernelParams );   //interpolation
		RBF( unsigned inDim, kernelType aKernel, Array<double> kernelParams, unsigned nKernel, double regParam, double maxErr, errMethod anErrMethod, unsigned maxKMeanIter, unsigned maxRegIter);	//supervised + ridge regression
		RBF( unsigned inDim, kernelType aKernel, Array<double> kernelParams, unsigned maximumBasis,  double maxErr);		//unsupervised , OLS-forward selection
		~RBF( );
		int evaluate( Individual& offspring );
		int evaluate( Population& offsprings );
		int evaluate( Array < double > inputData, Array < double >& outputData );
		double mse( Database Data );
		double mse( Array < double > Input, Array < double > Target );
		double mse( Population offsprings );
		int train( Array < double > InputData, Array < double > TargetData );
		int train( Database Data );

		void setDesignMatrix( Array<double>& inputData);
		void setWeightMatrix( Array<double>& targetData);
		Array<double> getDesignMatrix();
		Array<double> getWeightMatrix();
		Array<double> getOrthogonalMatrix();
		Array<double> getOrthogonalWeightMatrix();
		Array<double> getBaseCentres();


		static double getDistance(Array<double> input, unsigned i, Array<double> centre, unsigned j);
		static double gaussFunc(Array<double> input, unsigned i, Array<double> centre, unsigned j, Array<double> params);
		static double cubicFunc(Array<double> input, unsigned i, Array<double> centre, unsigned j, Array<double> params);
		static double linearFunc(Array<double> input, unsigned i, Array<double> centre, unsigned j, Array<double> params);
		static double thinPlateFunc(Array<double> input, unsigned i, Array<double> centre, unsigned j, Array<double> params);
		static double multiquadricFunc(Array<double> input, unsigned i, Array<double> centre, unsigned j, Array<double> params);
		static double inverseMultiquadricFunc(Array<double> input, unsigned i, Array<double> centre, unsigned j, Array<double> params);
		static double cauchyFunc(Array<double> input, unsigned i, Array<double> centre, unsigned j, Array<double> params);
		
		static double YPYFunc(Array<double> desMatrix, Array<double> output, double regParam, Array<double>& A); //Sum Squared Error
		static double GCVFunc(Array<double> desMatrix, Array<double> output, double regParam); // Generalized Cross-Validation
		static double UEVFunc(Array<double> desMatrix, Array<double> output, double regParam); // Unbiased Estimate Variance
		static double FPEFunc(Array<double> desMatrix, Array<double> output, double regParam); // Final Prediction Error
		static double BICFunc(Array<double> desMatrix, Array<double> output, double regParam); // Bayesian Information Criterion

		static double estRegParamGCV(double regParam, Array<double> weightMatrix, Array<double> designMatrix, Array<double>predicted);
		static double estRegParamUEV(double regParam, Array<double> weightMatrix, Array<double> designMatrix, Array<double>predicted);
		static double estRegParamFPE(double regParam, Array<double> weightMatrix, Array<double> designMatrix, Array<double>predicted);
		static double estRegParamBIC(double regParam, Array<double> weightMatrix, Array<double> designMatrix, Array<double>predicted);

		static void KMean(unsigned nCluster, unsigned maxIter, Array<double> input, Array<double>& centresFound, Array<unsigned>& clustMap, ofstream& fileStream);

		template <class Type>
		static void initArrayToZero( Array<Type>& arr ) {
		        unsigned i,j;
		        for(i = 0; i < arr.dim(0); i++) {
                		for(j=0; j<arr.dim(1); j++) {
                        		arr(i,j) = 0;
                		}
        		}
		}
	
};






