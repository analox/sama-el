#define BASE_DMODEL_ON_DF
//===========================================================================
/*!
 *  \file FFNet.cpp
 *
 *  \brief Offers the functions to create and to work with 
 *         a feed-forward network.
 *
 *  \author  C. Igel
 *  \date    2002
 *
 *  \par Copyright (c) 2002-2001:
 *      Institut f&uuml;r Neuroinformatik<BR>
 *      Ruhr-Universit&auml;t Bochum<BR>
 *      D-44780 Bochum, Germany<BR>
 *      Phone: +49-234-32-25558<BR>
 *      Fax:   +49-234-32-14209<BR>
 *      eMail: Shark-admin@neuroinformatik.ruhr-uni-bochum.de<BR>
 *      www:   http://www.neuroinformatik.ruhr-uni-bochum.de<BR>
 *
 *  \par Project:
 *      ReClaM
 *
 *  \par Language and Compiler:
 *      C++, egcs (Linux)
 *
 *  This file is part of ReClaM. This library is free software;
 *  you can redistribute it and/or modify it under the terms of the
 *  GNU General Public License as published by the Free Software
 *  Foundation; either version 2, or (at your option) any later version.
 *
 *  This library is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
 *  GNU General Public License for more details.
 *
 *  You should have received a copy of the GNU General Public License
 *  along with this library; if not, write to the Free Software
 *  Foundation, Inc., 675 Mass Ave, Cambridge, MA 02139, USA.
 *
 *
 */
//===========================================================================

#include "ReClaM/FFNet.h"

using namespace std;

//===========================================================================
/*!
 *  \brief Destructs a feed forward network object.
 *
 *  This destructor is used to replace the standard destructor,
 *  created by the compiler, because internal pointer
 *  variables must be deleted, too.
 *
 *  \warning none
 *  \bug     none
 *
 *  \author  C. Igel
 *  \date    2002
 *
 *  \par Changes
 *      none
 *
 *  \par Status
 *      none
 *
 */
FFNet::~FFNet() {
  if(z) delete [] z;
}

//===========================================================================
/*!
 *  \brief Activation function \f$g_{hidden}(x)\f$ of the hidden neurons.
 *
 *  The activation function is used for the propagation of the input
 *  through the network.
 *  Given a network with \f$M\f$ neurons, including \f$c\f$
 *  input neurons and \f$n\f$ output neurons, the sigmoid activation
 *  function for the hidden neuron with index
 *  \f$i \mbox{,\ } c \leq i < M - n\f$    is given as
 *
 *  \f$
 *      z_i = g_{hidden}(x) = \frac{1}{1 + \exp (-a)} 
 *  \f$
 *
 *  where \f$a\f$ as the propagated result of the input for
 *  the previous neurons is calculated as 
 *
 *  \f$
 *      a = \left( \sum_{j=0}^{j<i} w_{ij} z_j + \Theta_i \right)
 *  \f$
 *
 *  and \f$ \Theta_i \f$ denotes the bias term.
 *
 *      \param  a Input for the activation function, see above.
 *      \return \f$ z_i \f$.
 *
 *  \author  C. Igel
 *  \date    2002
 *
 *  \par Changes
 *      none
 *
 *  \par Status
 *      stable
 *
 *  \sa gOutput
 *
 */  
double FFNet::g(double a) {
  return 1 / (1 + exp(-a));
}


//===========================================================================
/*!
 *  \brief Computes the derivative of the activation function
 *         \f$g_{hidden}(x)\f$ for the hidden neurons.
 *
 *  The derivative function \f$g^*_{\mbox{hidden}}\f$ is defined
 *  as 
 * 
 *  \f$
 *      g^*_{hidden}(g_{hidden}(x)) = 
 *      \frac{\partial g_{hidden}(x)}{\partial x}
 *  \f$ 
 *      \param  ga The value of \f$g_{hidden}(x)\f$.
 *      \return The result of \f$g^*_{hidden}(g_{hidden}(x))\f$
 *
 *  \author  C. Igel
 *  \date    2002
 *
 *  \par Changes
 *      none
 *
 *  \par Status
 *      stable
 *
 *  \sa g
 *
 */  
double FFNet::dg(double ga) {
  return ga * (1 - ga);
}


//===========================================================================
/*!
 *  \brief Activation function \f$g_{output}(x)\f$ of the output neurons.
 *
 *  The activation function is used for the propagation of the input
 *  through the network.
 *  Given a network with \f$M\f$ neurons, including \f$c\f$
 *  input neurons and \f$n\f$ output neurons, the sigmoid activation
 *  function for the output neuron with index 
 *  \f$i \mbox{,\ } M - n \leq i < M \f$ 
 *  is given as
 *
 *  \f$
 *      z_i = g_{output}(x) = \frac{1}{1 + \exp (-a)} 
 *  \f$
 *
 *  where \f$a\f$ as the propagated result of the input for
 *  the previous neurons is calculated as 
 *
 *  \f$
 *      a = \left( \sum_{j=0}^{j<i} w_{ij} z_j + \Theta_i \right) 
 *  \f$
 *
 *  and \f$ \Theta_i \f$ denotes the bias term.
 *
 *  \param  a Input for the activation function, see above.
 *  \return \f$ z_i \f$.
 *
 *  \author  C. Igel
 *  \date    2002
 *
 *  \par Changes
 *      none
 *
 *  \par Status
 *      stable
 *
 *  \sa g
 *
 */  
double FFNet::gOutput(double a) {
  return 1 / (1 + exp(-a));
}

//===========================================================================
/*!
 *  \brief Computes the derivative of the activation function
 *         \f$g_{output}(x)\f$ for the output neurons.
 *
 *  The derivative function \f$g^*_{\mbox{output}}\f$ is defined
 *  as 
 * 
 *  \f$
 *      g^*_{output}(g_{output}(x)) = 
 *      \frac{\partial g_{output}(x)}{\partial x}
 *  \f$ 
 *
 *  \param  ga The value of \f$g_{output}(x)\f$.
 *  \return The result of \f$g^*_{output}(g_{output}(x))\f$
 *
 *  \author  C. Igel
 *  \date    2002
 *
 *  \par Changes
 *      none
 *
 *  \par Status
 *      stable
 *
 *  \sa gOutput, dg
 *
 */  
double FFNet::dgOutput(double ga) {
  return ga * (1 - ga);
}

//===========================================================================
/*!
 *  \brief Initializes the random number generator and the weights of the net.
 *
 *  After initializing the random number generator with a \em seed value,
 *  the weight values of the network are initialized with uniformly 
 *  distributed random numbers.
 *
 *      \param  seed Initialization value for random number generator, default
 *                   value is "42".
 *      \param  l Lower bound for weight random numbers, default value
 *                is "-0.5".
 *      \param  h Upper bound for weight random numbers, default value
 *                is "0.5".
 *      \return None.
 *
 *  \author  C. Igel
 *  \date    2002
 *
 *  \par Changes
 *      none
 *
 *  \par Status
 *      stable
 *
 */  
void FFNet::initWeights(long seed, double l, double h ) {
  Rng::seed(seed);
  for(unsigned i = inputDimension; i < numberOfNeurons; i++) {
    for(unsigned j = 0; j < i; j++) 
      if(connection(i, j)) weight(i, j) = Rng::uni(l, h); 
    if(connectionMatrix(i, bias)) weight(i, bias) = Rng::uni(l, h);
  }
  writeParameters();
}

//===========================================================================
/*!
 *  \brief Initializes the weights of the net.
 *
 *  The weight values of the network are initialized with uniformly 
 *  distributed random numbers.
 *
 *      \param  l Lower bound for weight random numbers, default value
 *                is "-0.5".
 *      \param  h Upper bound for weight random numbers, default value
 *                is "0.5".
 *      \return None.
 *
 *  \author  C. Igel
 *  \date    2002
 *
 *  \par Changes
 *      none
 *
 *  \par Status
 *      stable
 *
 *  \sa initWeights(long seed, double l, double h)
 *
 */  
inline 
void FFNet::initWeights(double l, double h) {
  for(unsigned i = inputDimension; i < numberOfNeurons; i++) {
    for(unsigned j = 0; j < i; j++) 
      if(connection(i, j)) weight(i, j) = Rng::uni(l, h); 
    if(connectionMatrix(i, bias)) weight(i, bias) = Rng::uni(l, h);
  }
  writeParameters();
}

//===========================================================================
/*!
 *  \brief Uses the input pattern(s) "in" for the Feed Forward Net
 *         model to produce the output vector(s) "output".
 *
 *  The given input pattern(s) in \em input are propagated forward
 *  through the net by using the #activate function for all
 *  patterns sequentially. After each propagation the result
 *  in the output neurons is stored in the output vector(s) \em output.
 *
 *  \param  input Input pattern(s) for the model.
 *  \param  output Output vector(s).
 *  \return None.
 *
 *  \author  C. Igel
 *  \date    2002
 *
 *  \par Changes
 *      none
 *
 *  \par Status
 *      stable
 *
 *  \sa activate
 *
 */  
inline
void FFNet::model(const Array<double> &input, Array<double> &output) {
  readParameters();
  if(input.ndim() == 1) {
    activate(input);
    for(unsigned i = firstOutputNeuron, c = 0 ; i < numberOfNeurons; i++, c++)
      output(c) = z[i];
  } else {
    for(unsigned pattern = 0; pattern < input.dim(0); pattern++) {
      activate(input[pattern]);
      for(unsigned i = firstOutputNeuron, c = 0 ; i < numberOfNeurons; i++, c++)
	output(pattern, c) = z[i];
    }
  }
}

//===========================================================================
/*!
 *  \brief Performs a backpropagation to calculate the derivatives of the 
 *         outputs with respect to the parameters of the network.
 *
 *  The values of the output neurons that were produced by an input pattern
 *  before are used for a backpropagation to calculate the partial derivatives 
 *  of the model function \f$f(x,w)\f$ with respect to the
 *  parameters 
 * 
 *  \f$
 *  \nabla f(w) = \left( \frac{\partial f}{\partial w_1}, \dots , 
 *                \frac{\partial f}{\partial w_n} \right)^T
 *  \f$
 *
 *  \return None.
 *
 *  \author  C. Igel
 *  \date    2002
 *
 *  \par Changes
 *      none
 *
 *  \par Status
 *      stable
 *
 */  
inline
void FFNet::backprop() {
  unsigned i, j;  // denote neurons
  unsigned c;     // denotes output
  unsigned pos;   // the actual position in the Array dmdw
  double sum;

  // 1. the 'local gradient' deltas for each
  //    output are calculated
  
  // output neurons
  for(j = firstOutputNeuron; j < numberOfNeurons; j++)
    for(c = 0; c < outputDimension; c++)
      if(c + firstOutputNeuron == j)  d(c, j) = dgOutput(z[j]);
      else d(c, j) = 0;
  
  // hidden units
  for(j = firstOutputNeuron - 1; j >= inputDimension; j--) {
    for(c = 0; c < outputDimension; c++) {
      sum = 0;
      for(i = j + 1; i < numberOfNeurons; i++) {
	if(connection(i, j)) sum += weight(i, j) * d(c, i);
      }
      d(c, j) = dg(z[j]) * sum;
    }
  }
  
  // 2. the gradient for each output is calculated from
  //    delta
  pos = 0;
  for(i = inputDimension; i < numberOfNeurons; i++) {
    for(j = 0; j < i; j++) {
      if(connection(i, j)) {
	for(c = 0; c < outputDimension; c++) { 
	  dmdw(c, pos) += d(c, i) * z[j];
	}
	pos++;
      }
    }
    if(connection(i, bias)) {
      for(c = 0; c < outputDimension; c++) { // bias 
	dmdw(c, pos) += d(c, i);
      }
      pos++;
    }
  }
}

//===========================================================================
/*!
 *  \brief Reads in one input pattern for the Feed Forward Net
 *         model and calculates the derivatives of the resulting network 
 *         output with respect to the weights. Furthermore, the
 *         network output is given back.
 *
 *  Equal to method #dmodel(const Array<double> &input), but here
 *  the output of the network is not only used for the calculation,
 *  but also stored in \em output.
 *
 *  \param  input Input pattern for the model.
 *  \param  output Output for the input pattern.
 *  \return None.
 *
 *  \author  C. Igel
 *  \date    2002
 *
 *  \par Changes
 *      none
 *
 *  \par Status
 *      stable
 *
 */  
void FFNet::dmodel(const Array<double> &input, Array<double> &output) {
  unsigned i;  // denotes neurons
  unsigned c;  // denotes output

  dmdw =0.;

#ifdef BASE_DMODEL_ON_DF
  Array<double> coeffs;
  coeffs.resize(outputDimension, false);
  Array<double> gradbuf;
  gradbuf.resize(w.nelem(), false);
#endif

  readParameters();

  if(input.ndim() == 1) {
    // calculate output
    activate(input);
    for(i = firstOutputNeuron, c = 0 ; i < numberOfNeurons; i++, c++)
      output(c) = z[i];
#ifdef BASE_DMODEL_ON_DF
    for(i = 0 ; i < outputDimension; i++) {
      coeffs = 0.; coeffs(i) = 1.;
      gradbuf = 0.;
      df(coeffs, gradbuf);
      dmdw[i] -= gradbuf;
    }
#else
    // backpropagation step
    backprop();
#endif
  } else {
    for(unsigned pattern = 0; pattern < input.dim(0); pattern++) {
      activate(input[pattern]);
      for(i = firstOutputNeuron, c = 0 ; i < numberOfNeurons; i++, c++)
	output(pattern, c) = z[i];
#ifdef BASE_DMODEL_ON_DF
      for(i = 0 ; i < outputDimension; i++) {
	coeffs = 0.; coeffs(i) = 1.;
	gradbuf = 0.;
	df(coeffs, gradbuf);
	dmdw[i] -= gradbuf;
      }
#else
      backprop();
#endif
    }
  }
};

inline void FFNet::df(const Array<double> &coeffs, Array<double> &dfdw) {
  unsigned i, j;  // denote neurons
  unsigned c;     // denotes output
  unsigned pos;   // the actual position in the Array dedw
  double sum;

  for(i = firstOutputNeuron, c = 0 ; i < numberOfNeurons; i++, c++) {
    delta(i) =  coeffs(c) * dgOutput(z[i]);
  }
  
  // hidden units
  for(j = firstOutputNeuron - 1; j >= inputDimension; j--) {
    sum = 0;
    for(i = j + 1; i < numberOfNeurons; i++) {
      if(connection(i, j)) sum += weight(i, j) * delta(i);
    }
    delta(j) = dg(z[j]) * sum;
  }
  
  // calculate error gradient
  pos = 0;
  for(i = inputDimension; i < numberOfNeurons; i++) {
    for(j = 0; j < i; j++) {
      if(connection(i, j)) {
	dfdw(pos) -= delta(i) * z[j];
	pos++;
      }
    }
    if(connection(i, bias)) {  // bias
      dfdw(pos) -= delta(i); 
      pos++;
    }
  }
};

//===========================================================================
/*!
 *  \brief Reads in one input pattern for the Feed Forward Net
 *         model and calculates the derivatives of the resulting network 
 *         output with respect to the weights.
 *
 *  The single input pattern \em in is used to #activate the neurons
 *  of the network. The results, given in the output neurons
 *  are then used to calculate the derivatives of the output
 *  with respect to the weights of the network.
 *
 *  \param  input The single input pattern. If more than one
 *                pattern is given, the method exits with
 *                failure.
 *  \return None.
 *
 *  \author  C. Igel
 *  \date    2002
 *
 *  \par Changes
 *      none
 *
 *  \par Status
 *      stable
 *
 *  \sa backprop
 *
 */  
inline 
void FFNet::dmodel  (const Array<double> &input) {
  readParameters();

  if(input.ndim() == 1) {
    activate(input);
    // backpropagation step
    backprop();
  } else {
    cerr << "the derivative of the model with respect to more" 
	 << " than one input is not defined" << endl;
    exit(EXIT_FAILURE);
  }
};


//===========================================================================
/*!
 *  \brief Reserves memory for all internal net data structures. 
 *
 *  For the internal administration of a network several dynamic
 *  data structures are used. Based on the size of the #connectionMatrix
 *  and the values of ModelInterface::inputDimension and
 *  ModelInterface::outputDimension, the memory for all other data
 *  is reserved.
 *
 *  \return None.
 *
 *  \author  C. Igel
 *  \date    2002
 *
 *  \par Changes
 *      2002-01-08, ci <br>
 *      Some commands added to remove memory leackage,
 *      that lead to crashes.
 *
 *  \par Status
 *      stable
 *
 */  
void FFNet::resize() {
  numberOfNeurons   = connectionMatrix.dim(0);
  bias              = numberOfNeurons;
  firstOutputNeuron = numberOfNeurons - outputDimension;

  unsigned numberOfWeights = 0; //numberOfNeurons - inputDimension;
  for(unsigned i = 0; i < connectionMatrix.dim(0); i++) {
    for(unsigned j = 0; j < i; j++) 
      if(connectionMatrix(i, j)) numberOfWeights++;
    if(connectionMatrix(i, bias)) numberOfWeights++;
  }

  w.resize (numberOfWeights);
  d.resize (outputDimension, numberOfNeurons);
  dmdw.resize(outputDimension, numberOfWeights);
  dedw.resize(numberOfWeights);
  delta.resize(numberOfNeurons, false);

  if (z) delete [] z;
  if (numberOfNeurons) z = new double[numberOfNeurons];
  else z = NULL;
}


//===========================================================================
/*!
 *  Constructor no. 1
 *
 *  \brief Creates an empty feed-forward network with "in" 
 *         input neurons and "out" output neurons.
 *
 *  Only the input and output dimensions are set, but the network
 *  will contain no neurons.
 *
 *      \param in Dimension of the input (no. of input neurons).
 *      \param out Dimension of the output (no. of output neurons).
 *      \return None.
 *
 *  \author  C. Igel
 *  \date    2002
 *
 *  \par Changes
 *      none
 *
 *  \par Status
 *      stable
 *
 */  
FFNet::FFNet(const unsigned in, const unsigned out) {
  inputDimension   = in;
  outputDimension  = out;

  numberOfNeurons   = 0;
  bias              = 0;
  firstOutputNeuron = 0;

  z = NULL;
}


//===========================================================================
/*!
 *  Constructor no. 2
 *
 *  \brief Creates a feed-forward network with "in" input neurons and
 *         "out" output neurons. Additionally, the array "cmat" determines
 *         the topology (i.e., number of neurons and their connections).
 *         
 *
 *  A network with the given connections will be created, memory for
 *  the #weightMatrix reserved, but the weights for all connections
 *  will be set to zero.
 *
 *      \param in Dimension of the input (no. of input neurons).
 *      \param out Dimension of the output (no. of output neurons).
 *      \param cmat The #connectionMatrix.
 *      \return None.
 *
 *  \author  C. Igel
 *  \date    2002
 *
 *  \par Changes
 *      2002-01-09, ci <br>
 *      Memory leackage removed.
 *
 *  \par Status
 *      stable
 *
 */  
FFNet::FFNet(const unsigned in, unsigned out, 
	     const Array<int>& cmat) {
  if(cmat.ndim()!=2) {
    cerr << "connection matrix has not dimension two" << endl;
    exit(EXIT_FAILURE);
  }
  if(cmat.dim(0) != (cmat.dim(1) - 1)) {
    cerr << "connection matrix has to be (N + 1) x N" << endl;
    exit(EXIT_FAILURE);
  }

  inputDimension   = in;
  outputDimension  = out;
  connectionMatrix = cmat;

  weightMatrix.resize(connectionMatrix); 
  weightMatrix     = 0;

  z = NULL;

  resize();
  writeParameters();
}


//===========================================================================
/*!
 *  Constructor no. 3
 *
 *  \brief Creates a feed-forward network with "in" input neurons and
 *         "out" output neurons. Additionally, the arrays "cmat" and
 *         "wmat" determine the topology (i.e., number of neurons and their
 *         connections) as well as the connection weights.
 *
 *  A network with the given connections and weights will be created.
 *
 *      \param in Dimension of the input (no. of input neurons).
 *      \param out Dimension of the output (no. of output neurons).
 *      \param cmat The #connectionMatrix.
 *      \param wmat The #weightMatrix.
 *      \return None.
 *
 *  \author  C. Igel
 *  \date    2002
 *
 *  \par Changes
 *      2002-01-09, ci <br>
 *      Memory leackage removed.
 *
 *  \par Status
 *      stable
 *
 */  
FFNet::FFNet(const unsigned in, unsigned out, 
	     const Array<int>& cmat, const Array<double>& wmat) {
  if(cmat.ndim()!=2) {
    cerr << "connection matrix has not dimension two" << endl;
    exit(EXIT_FAILURE);
  }
  if(cmat.dim(0) != (cmat.dim(1) - 1)) {
    cerr << "connection matrix has to be (N + 1) x N" << endl;
    exit(EXIT_FAILURE);
  }
  if(wmat.ndim()!=2) {
    cerr << "weight matrix has not dimension two" << endl;
    exit(EXIT_FAILURE);
  }
  if(wmat.dim(0) != (wmat.dim(1) - 1)) {
    cerr << "weight matrix has to be (N + 1) x N" << endl;
    exit(EXIT_FAILURE);
  }

  inputDimension   = in;
  outputDimension  = out;
  connectionMatrix = cmat;
  weightMatrix     = wmat; 

  z = NULL;
  resize();
  writeParameters();
}


//===========================================================================
/*!
 *  Constructor no. 4
 *
 *  \brief Creates a feed-forward network by reading the necessary 
 *         information from a file named "filename".
 *
 *  A file is used to create a new feed-forward network. This
 *  file must have the following content:
 *  The first line of the file must contain two numbers that specify 
 *  the number of input and the number of output neurons of the
 *  network, respectively.<BR>
 *  This line is followed by the values for the #connectionMatrix.<BR>
 *  The third and last part are the values for the #weightMatrix. 
 *
 *  \param filename Name of the file that contains the information
 *                  for the creation of the network. If the file
 *                  doesn't exist, the method will exit with
 *                  failure.
 *  \return None.
 *
 *  \par Example
 *  <BR>
 *  1 1<BR>
 *  <BR> 
 *  0 0 0<BR>
 *  1 0 0<BR>
 *  0 1 0<BR>
 *  <BR>  
 *  0 0 0 2<BR>
 *  3 0 0 2<BR>
 *  0 3 0 2<BR>
 *  <BR>
 *
 *  A file with the content shown above will create a network
 *  with 1 input and 1 output neuron.<BR>
 *  A connection exists from the input neuron to the single
 *  hidden neuron of the network and from the hidden neuron
 *  to the output neuron. Each of the two connections has
 *  a weight of "3".<BR>
 *  The connection of each neuron to the bias value has a weight of "2".
 *
 *  \author  C. Igel
 *  \date    2002
 *
 *  \par Changes
 *      none
 *
 *  \par Status
 *      stable
 *
 */  
FFNet::FFNet(const string &filename) {
  ifstream input(filename.c_str());
  if(!input) {
    cerr << "cannot open net file " << filename << endl;
    exit(EXIT_FAILURE);
  }
  z = NULL;
  input >> *this;
  input.close();
}


//===========================================================================
/*!
 *  \brief Based on a given #connectionMatrix and a #weightMatrix the
 *         structure of the current network is modified.
 *
 *  This function will use the information in the
 *  given #connectionMatrix \em cmat and the #weightMatrix \em wmat to 
 *  modify the network.
 *
 *      \param cmat The #connectionMatrix that provides the
 *                  basis information for the modification
 *                  of the network.
 *      \param wmat The #weightMatrix that provides the
 *                  basis information for the modification
 *                  of the network.
 *      \return     None.
 *
 *  \author  C. Igel
 *  \date    2002
 *
 *  \par Changes
 *      none
 *
 *  \par Status
 *      stable
 *
 */  
void FFNet::setStructure(const Array<int>& cmat, const Array<double>& wmat) {
  if((cmat.ndim()!=2) || 
     (wmat.ndim()!=2)) {
    cerr << "connection or weight matrix has not dimension two" << endl;
    exit(EXIT_FAILURE);
  }

  if((cmat.dim(0) != (cmat.dim(1) - 1)) || 
     (cmat.dim(0) != wmat.dim(0)) ||
     (cmat.dim(1) != wmat.dim(1))) {
    cerr << "connection and weight matrices have to be (N + 1) x N" << endl;
    exit(EXIT_FAILURE);
  }

  connectionMatrix = cmat;
  weightMatrix     = wmat;

  resize();
  writeParameters();
}    

//===========================================================================
/*!
 *  \brief Based on a given #connectionMatrix the structure of the
 *         current network is modified.
 *
 *  This function will use the information in the
 *  given #connectionMatrix \em cmat to modify the network.
 *  The weights for all neurons of the new #connectionMatrix
 *  are initialized by zero.
 *
 *      \param cmat The #connectionMatrix that provides the
 *                  basis information for the creation
 *                  of the network.
 *      \return     None.
 *
 *  \author  C. Igel
 *  \date    2002
 *
 *  \par Changes
 *      none
 *
 *  \par Status
 *      stable
 *
 */  
void FFNet::setStructure(const Array<int>& cmat) {
  if(cmat.ndim()!=2) {
    cerr << "connection matrix has not dimension two" << endl;
    exit(EXIT_FAILURE);
  }
  if(cmat.dim(0) != (cmat.dim(1) - 1)) {
    cerr << "connection matrix has to be (N + 1) x N" << endl;
    exit(EXIT_FAILURE);
  }

  connectionMatrix = cmat;
  weightMatrix.resize(cmat.dim(0), cmat.dim(1), false);

  resize();
  w = 0;

  writeParameters();
}    


//===========================================================================
/*!
 *  \brief Based on a given #weightMatrix the structure of the current
 *         network is modified.
 *
 *  This function will use the values of ModelInterface::inputDimension,
 *  ModelInterface::outputDimension and the information in the
 *  given #weightMatrix \em wmat to change the structure of the network.
 *  An existing connection matrix (that is compatible to the 
 *  the new #weightMatrix \em wmat) can be adopted or a new one
 *  will be created with the data of the #weightMatrix.
 *
 *      \param wmat     The #weightMatrix that provides the
 *                      basis information for the modification
 *                      of the network.
 *      \param preserve If set to "true", the existing
 *                      #connectionMatrix will be kept.
 *                      In case of inconsistency between
 *                      the #weightMatrix and the #connectionMatrix
 *                      the function will exit with failure.
 *      \return         None.
 *
 *  \author  C. Igel
 *  \date    2002
 *
 *  \par Changes
 *      2002-01-08, ci <br>
 *      Added some commands to remove memory leackage,
 *      that lead to crashes, too.
 *
 *  \par Status
 *      stable
 *
 */  
void FFNet::setStructure(const Array<double>& wmat, bool preserve) {
  if(wmat.ndim()!=2) {
    cerr << "weight matrix has not dimension two" << endl;
    exit(EXIT_FAILURE);
  }
  if(wmat.dim(0) != (wmat.dim(1) - 1)) {
    cerr << "weight matrix has to be (N + 1) x N" << endl;
    exit(EXIT_FAILURE);
  }

  weightMatrix = wmat;
  connectionMatrix.resize(wmat.dim(0), wmat.dim(1), false);

  numberOfNeurons   = connectionMatrix.dim(0);
  bias              = numberOfNeurons;
  firstOutputNeuron = numberOfNeurons - outputDimension;

  unsigned numberOfWeights = 0;
  for(unsigned i = 0; i < connectionMatrix.dim(0); i++) {
    for(unsigned j = 0; j < i; j++) {
      if(preserve) {
	if(weightMatrix(i, j)  != 0.) {
	  if(connectionMatrix(i, j) == 0) {
	    cerr << "weight specified for non-existing connection" << endl;
	    exit(EXIT_FAILURE);
	  }
	  numberOfWeights++;
	} else if(connectionMatrix(i, j) != 0) numberOfWeights++;
	if(weightMatrix(i, bias)  != 0.) {
	  if(connectionMatrix(i, bias) == 0) {
	    cerr << "weight specified for non-existing bias" << endl;
	    exit(EXIT_FAILURE);
	  }
	  numberOfWeights++;
	} else if(connectionMatrix(i, bias) != 0) numberOfWeights++;
      } else {
	if(weightMatrix(i, j)  != 0.) {
	  connectionMatrix(i, j) = 1;
	  numberOfWeights++;
	} else {
	  connectionMatrix(i, j) = 0;	  
	}
	if(weightMatrix(i, bias)  != 0.) {
	  connectionMatrix(i, bias) = 1;
	  numberOfWeights++;
	} else {
	  connectionMatrix(i, bias) = 0;	  
	}
      }
    }
  }
  
  resize ();
  writeParameters();
}    

//===========================================================================
/*!
 *  \brief Reads the parameters (weights) from the 
 *         #weightMatrix and writes them to the parameter vector
 *         ModelInterface::w.
 *
 *  In feed forward networks connections and the corresponding weights
 *  can exist only from neurons with no. \f$j\f$ to neurons with no. \f$i\f$ 
 *  with \f$j < i\f$.
 *  So to save space, the parameters, i.e. the values for the connection
 *  weights are internally stored in a vector.
 *  But a matrix offers a more concise view so this function reads in the 
 *  original (extended) #weightMatrix and stores the values in the
 *  parameter vector ModelInterface::w.
 *
 *  \return None.
 *
 *  \author  C. Igel
 *  \date    2002
 *
 *  \par Changes
 *      none
 *
 *  \par Status
 *      stable
 *
 */  
void FFNet::writeParameters() {
  for(unsigned k = 0, i = 0; i < connectionMatrix.dim(0); i++) {
    for(unsigned j = 0; j < i; j++) {
      if(connectionMatrix(i, j)) {
	w(k) = weight(i, j);
	k++;
      }
    }
    if(i >= inputDimension) {
      if(connectionMatrix(i, bias)) {
	w(k) = weight(i, bias);
	k++;
      }
    }
  }
}  

//===========================================================================
/*!
 *  \brief Reads the parameters (weights) from the parameter
 *         vector ModelInterace::w and stores them in the #weightMatrix.
 *
 *  In feed forward networks connections and the corresponding weights
 *  can exist only from neurons with no. \f$j\f$ to neurons with no. 
 *  \f$i\f$ with \f$j < i\f$.
 *  So to save space, the parameters, i.e. the values for the connection
 *  weights are internally stored in a vector.
 *  But a matrix offers a more concise view so this function reads in the 
 *  parameter vector ModelInterface::w and writes the values to the  
 *  #weightMatrix.
 *
 *  \return None.
 *
 *  \author  C. Igel
 *  \date    2002
 *
 *  \par Changes
 *      none
 *
 *  \par Status
 *      stable
 *
 */  
void FFNet::readParameters() {
  for(unsigned k = 0, i = 0; i < connectionMatrix.dim(0); i++) {
    for(unsigned j = 0; j < i; j++) {
      if(connectionMatrix(i, j)) {
	weight(i, j) = w(k);
	k++;
      }
    }
    if(i >= inputDimension) {
      if(connectionMatrix(i, bias)) {
	weight(i, bias) = w(k);
	k++;
      }
    }
  }
}  


//===========================================================================
/*!
 *  \brief Returns the current #weightMatrix.
 *
 *  The weight values are read from the parameter vector
 *  ModelInterface::w and stored in the #weightMatrix,
 *  then the matrix is returned.
 *
 *  \return The #weightMatrix.
 *
 *  \author  C. Igel
 *  \date    2002
 *
 *  \par Changes
 *      none
 *
 *  \par Status
 *      stable
 *
 */  
const Array< double > &FFNet::getWeights() {
  readParameters();
  return weightMatrix;
}

//===========================================================================
/*!
 *  \brief Returns the current #connectionMatrix.
 *
 *  \return The #connectionMatrix.
 *
 *  \author  C. Igel
 *  \date    2002
 *
 *  \par Changes
 *      none
 *
 *  \par Status
 *      stable
 *
 */  
const Array< int > &FFNet::getConnections() {
  return connectionMatrix;
}

//===========================================================================
/*!
 *  \brief Writes the data of a network to an output stream.
 *
 *  The ModelInterface::inputDimension, ModelInterface::outputDimension,
 *  the #connectionMatrix and the #weightMatrix are written to an output 
 *  stream.
 *  For the syntax of a possible output, see the example section of
 *  #FFNet(const std::string &filename).
 *
 *      \param os The output stream to where the data is written.
 *      \param net The net that will be written to the output stream.
 *      \return Reference to the output stream.
 *
 *  \author  C. Igel
 *  \date    2002
 *
 *  \par Changes
 *      none
 *
 *  \par Status
 *      stable
 *
 */  
ostream& operator<< (ostream& os, FFNet &net)
{
  os << net.inputDimension << " " << net.outputDimension << "\n";
  writeArray(net.connectionMatrix, os,"","\n",' ');
  writeArray(net.weightMatrix, os,"","\n",' ');
  return os;
}


//===========================================================================
/*!
 *  \brief Reads the data for the creation of a network from an input stream.
 *
 *  The ModelInterface::inputDimension, ModelInterface::outputDimension, 
 *  the #connectionMatrix and the
 *  #weightMatrix are read from the input stream,
 *  memory for the internal net data structures is reserved
 *  and the weight and bias values are written to the 
 *  parameter vector ModelInterface::w.
 *  For the syntax of a possible input, see the example section of
 *  #FFNet(const std::string &filename).
 *
 *      \param is The input stream, from which the data is read.
 *      \param net The net object, to which the data will be copied.
 *      \return Reference to the input stream.
 *
 *  \author  C. Igel
 *  \date    2002
 *
 *  \par Changes
 *      none
 *
 *  \par Status
 *      stable
 *
 */  
istream& operator>> (istream& is, FFNet &net)
{
  unsigned i, j, n, pos;
  double v, dn;
  vector<double> l;

  is >> net.inputDimension;
  is >> net.outputDimension;

  is >> v;
  while(is) {
    l.push_back(v);
    is >> v;
  }

  if(!l.size()) {
    cerr << "no matrices given in net file";
    exit(EXIT_FAILURE);
  }

  dn = (sqrt(1. + l.size() * 2.) - 1.) / 2.;
  n = unsigned(dn);
  if(double(n) == dn) { // new style with bias
    net.connectionMatrix.resize(n, n + 1);
    net.weightMatrix.resize(n, n + 1);

    pos = 0;
    for(i = 0; i < n; i++) {
      for(j = 0; j < n + 1; j++) {
	net.connectionMatrix(i, j) = unsigned(l[pos]);
	pos++;
      }
    }

    for(i = 0; i < n; i++) {
      for(j = 0; j < n + 1; j++) {
	net.weightMatrix(i, j) = l[pos];
	pos++;
      }
    }
  } else { // old style without bias
    dn = (sqrt(1. + l.size() * 8.) - 1.) / 4.;
    n =  unsigned(dn);
    if(double(n) == dn) {
      net.connectionMatrix.resize(n, n + 1);
      net.weightMatrix.resize(n, n + 1);

      pos = 0;
      for(i = 0; i < n; i++) {
	for(j = 0; j < n; j++) {
	  net.connectionMatrix(i, j) = unsigned(l[pos]);
	  pos++;
	}
	// handle bias
	if(i >= net.inputDimension) net.connectionMatrix(i, n) = 1; 
	else net.connectionMatrix(i, n) = 0;
      }
      for(i = 0; i < n; i++) {
	for(j = 0; j < n + 1; j++) {
	  net.weightMatrix(i, j) = l[pos];
	  pos++;
	}
      }
    } else { // no style
      cerr << "cannot parse network configuration file" << endl;
      exit(EXIT_FAILURE);
    }
  }

  net.resize();
  net.writeParameters();

  return is;
}



