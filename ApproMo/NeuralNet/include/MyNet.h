#ifndef __MyNet__
#define __MyNet__

#include <ReClaM/MSEFFNet.h>
#include <ReClaM/Rprop.h>
#include <ReClaM/MeanSquaredError.h>
#include <ReClaM/ErrorMeasures.h>

using namespace std;
//! MyNet generates an initial neural network model in different ways.
/*! This class is based on classes defined in ReClaM of the Shark library.
  The following ReClaM functions are used in \c metaModel.h and \c UpdatemetaModel.h .
  -# \c MyNet.model(inputNN,outputNN)
  This function returns the output of a neural
  network for a given input vector \c inputNN.
  The output value
  is stored in the variable \c outputNN. 
  -#  \c metaNet.initRprop(delta0),
  \c metaNet.rprop(inputNN,targetNN,np,nm,dmax,dmin); \n
  There exist several methods for training neural networks. 
  The training algorithm used in this software is the Rprop, which is
  a robust and fast training algorithm.  \n The parameters of the
  training algorithms are:
  - delta0: the initial learning rate
  - np: the increase factor if the learning rate should increase
  - nm: the decrease factor if the learning rate should decrease
  - dmax: the maximal increase factor allowed
  - dmin: the minimal increase factor allowed
  - inpuNN: the input vector  of the training data
  - targetNN: the output vector of the training data 
  -#  \c metaNet.inputDimension(),metaNet.outputDimension
  These two functions return the input dimension and output dimension
  of a neural network if it has been defined by a constructor.
  -#  \c metaNet.getConnections();
  This function returns the struture of the neural network which is
  specified by a connection matrix.
  -#  \c metaNet.getWeights()
  This function returns the weight for each connection in the neural
  network.
*/

class MyNet : public MSEFFNet,       // feed-forward network
  public IRpropPlus,     // resilent backpropagation
  public ErrorMeasures   // winner-takes-all & 
{
 public:
  //! Initialize a neural network model by reading in a file
  /*! Both the structure and parameters of the neural network
    are specified in the file. Information include the number of
    inputs, the number of hidden nodes and the number of outputs.
    Besides, the connection between the nodes are also determined
    by a matrix. If two nodes are connected, the corresponding
    element in the matrix is 1, otherwise, it is 0. This is the
    initialization mode adopted in the current software.*/
  MyNet      (const string &filename) :
    MSEFFNet(filename) {  }

  MyNet(const unsigned in = 0, const unsigned out =0):
    MSEFFNet(in, out) {}
 


  //! Initialize a neural network model by specifying the structure and weights separately
  /*! This function initializes an neural network model by giving the number
    of inputs, the number of outputs, a matrix that determines the connections
    between the neurons and a matrix specifying the weight of each connection.
  */
  MyNet(const unsigned in, const unsigned out,
        const Array<int> &cmat, const Array<double> &wmat) :
    MSEFFNet(in, out, cmat, wmat) {}


#ifdef _WINDOWS
  double error(const Array<double> &input, const Array<double> &target) {
    return MSEFFNet::error(input,target);
  }
  double derror(const Array<double> &input, const Array<double> &target,
		bool returnError = true) {
    return MSEFFNet::derror(input,target,returnError);
  }
  
  void dmodel(const Array<double> &input) {
    MSEFFNet::dmodel(input);
  }
  
  void dmodel(const Array<double> &input, Array<double> &output) {
    MSEFFNet::dmodel(input,output);
  }
  
  void   model         (const Array<double> &input, Array<double> &output) {
    FFNet::model(input,output);
    
  }
  
  void df(const Array<double> &coeffs, Array<double> &dfdw) {
    
    unsigned i, j;  // denote neurons
    unsigned c;     // denotes output
    unsigned pos;   // the actual position in the Array dedw
    double sum;
    
    for(i = firstOutputNeuron, c = 0 ; i < numberOfNeurons; i++, c++) {
      FFNet::delta(i) =  coeffs(c) * dgOutput(z[i]);
    }
    
    // hidden units
    for(j = firstOutputNeuron - 1; j >= inputDimension; j--) {
      sum = 0;
      for(i = j + 1; i < numberOfNeurons; i++) {
	if(connection(i, j)) sum += weight(i, j) * FFNet::delta(i);
      }
      FFNet::delta(j) = dg(z[j]) * sum;
    }
    
    // calculate error gradient
    pos = 0;
    for(i = inputDimension; i < numberOfNeurons; i++) {
      for(j = 0; j < i; j++) {
	if(connection(i, j)) {
	  dfdw(pos) -= FFNet::delta(i) * z[j];
	  pos++;
	}
      }
      if(connection(i, bias)) {  // bias
	dfdw(pos) -= FFNet::delta(i); 
	pos++;
      }
    }
  };
  
#endif
  
  
  
  
  //! Alternative activation functions for hidden nodes
  double g       (double a)  { return a / (1 + fabs(a)); }
  double dg      (double ga) { return  (1 - sgn(ga) * ga) * (1 - sgn(ga) * ga); };

  //! Linear output nodes
  double gOutput       (double a)  { return a; }
  double dgOutput      (double ga) { return 1; };
};
#endif
