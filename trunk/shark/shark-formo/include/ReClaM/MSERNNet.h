//===========================================================================
/*!
 *  \file MSERNNet.h
 *
 *  \brief Offers the functions to create and to work with a
 *         recurrent neural network.
 *         The network is combined with the mean squared error
 *         measure. This combination is created due to computational
 *         efficiency. 
 *         
 *  \author  M. Toussaint
 *  \date    2000
 *
 *  \par Copyright (c) 1999-2001:
 *      Institut f&uuml;r Neuroinformatik<BR>
 *      Ruhr-Universit&auml;t Bochum<BR>
 *      D-44780 Bochum, Germany<BR>
 *      Phone: +49-234-32-27974<BR>
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
 *  \par File and Revision:
 *      $RCSfile: MSERNNet.h,v $<BR>
 *      $Revision: 2.3 $<BR>
 *      $Date: 2005/08/03 12:04:30 $
 *
 *  \par Changes:
 *      $Log: MSERNNet.h,v $
 *      Revision 2.3  2005/08/03 12:04:30  christian_igel
 *      changes limits/values, MAXDOUBLE etc.
 *
 *      Revision 2.2  2004/06/16 15:22:20  saviapbe
 *
 *      Some bugs in the doxygen documentation were removed.
 *
 *      Revision 2.1  2004/05/27 15:13:58  saviapbe
 *      The documentation was changed.
 *
 *      Revision 2.0  2003/11/28 16:23:22  shark-admin
 *      Revision tag reset to revision tag 2.x
 *
 *      Revision 1.1.1.1  2003/11/24 13:37:06  shark-admin
 *      INI Administration
 *
 *      Revision 1.3  2002/05/16 13:26:03  rudi
 *      doxygen commands added/modified.
 *
 *      Revision 1.2  2002/02/06 14:47:38  rudi
 *      Doxygen comments added.
 *
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

#include <cstdio>
#include <fstream>
#include <iomanip>
#include "Array/ArrayTable.h"
#include "Array/ArrayOp.h"
#include "Rng/GlobalRng.h"
#include "ReClaM/ModelInterface.h"

#ifndef __SOALRIS__
#include <limits>
#else
#include <values.h>
#endif

#ifndef RNNET_H
#define RNNET_H

//===========================================================================
/*!
 *  \brief A recurrent neural network regression model that learns
 *         with Back Propagation Through Time and assumes the MSE
 *         error functional.
 *         
 *  This class defines a recurrent neural network regression
 *  model. The input and output arrays should, as usual, be
 *  2-dimensional and have the dimensionality
 *  (batch-length,number-of-neurons). The only difference is that
 *  here, the batch corresponds to a time series instead of
 *  independent data points. The gradient is calculated via
 *  BackPropTroughTime (BPTT). This model implies the MSE as the error
 *  functional.
 *
 *  The class can handle arbitrary network architectures. Feed-forward
 *  connections (with time delay `zero') as well as connections of any
 *  time delay can be realized --- see setStructure() for details of how
 *  to define the structure.
 *
 *  Things to note!:
 *
 *  (1) All neurons are sigmoidal! So please transform inputs and
 *  outputs to take this into account. The reason is that internally,
 *  no differences between input, hidden, or output neurons is made.
 *
 *  (2) The initialization of the state history is an important issue
 *  for dynamic systems like an RNN. Please read the documentation of
 *  setWarmUpLentgh().
 *
 *  (3) Online learning can, in principle, be realized by having a
 *  batch-length of zero (the input/output arrays still have to be
 *  2-dimensional, i.e., of dimension (1,number-of-neurons)) and by
 *  setting the WarmUpLength to zero. Note though that this is quite
 *  inefficient compared to batch-learning (because the BPTT does not
 *  save you any calculation time) and that cross-talk will make
 *  learning _very_ difficult. Maybe, even when online learning, you
 *  should still use a reasonable batch size (say 100), feed these
 *  batches sequentially, and set the WarmUpLength to zero.
 *
 *  Please also refer to the
 *
 *  \example simpleRNNet.cpp
 *
 *  (Internals: the only core functions are `processTimeSeries' and
 *  `calcGradBPTT' in the cpp-file. Almost all the rest it `utility'
 *  stuff.)
 *
 *  \author  M. Toussaint
 *  \date    2000
 *
 *  \par Changes:
 *      none
 *
 *  \par Status:
 *      stable
 */
class MSERNNet:virtual public ModelInterface{


// The following command is used for the
// doxygen documentation program. doxygen is forced to include the
// source code of the named example program into the documentation,
// that then can be linked directly.
// This "dirty trick" is necessary, because to put "\example [filename]"
// directly into method documentations has some unaesthetic side effects.
//
/*! \example simpleRNNet.cpp */


public:


  //! Creates an empty MSE recurrent neural network.
  MSERNNet();

  //! Creates an empty MSE recurrent neural network with a certain
  //! topology given by the connection matrix "con".
  MSERNNet(Array<int>);

  //! Creates an empty MSE recurrent neural network with a certain
  //! topology given by the connection matrix "con" and a
  //! set of weights given by weight matrix "wei".
  MSERNNet(Array<int>,Array<double>);

  //! Creates a recurrent neural network by reading the necessary 
  //! information from a file named "filename".
  MSERNNet(std::string);

  //! Sets the length of the warmup sequence.
  void includeWarmUp(unsigned =0);

  //! Propagates input patterns through the network to modify 
  //! internal states.
  void model(const Array<double> & input);

  //! Propagates input patterns through the network to modify internal 
  //! states and to determine the network output.
  void model(const Array<double> & input,Array<double> &output);

  //! Evaluates the error percentage of the network output compared 
  //! to given target data.
  double errorPercentage(const Array<double> &input,const Array<double> &target);

  //! Evaluates the mean squared error of the network output compared 
  //! to given target data.
  double error(const Array<double> &input,const Array<double> &target);

  //! Evaluates the derivative of the mean squared error.
  double derror(const Array<double> &input,const Array<double> &target, 
                bool returnError = true);

  //! Evaluates the mean squared error of the network output compared 
  //! to given target data.
  double meanSquaredError(const Array<double> &input,const Array<double> &target);

  //! Returns the connection Matrix.  
  Array<int> &getConnections();

  //! Returns the weight matrix.
  Array<double> &getWeights();

  //! Initializes the weights of the net.
  void initWeights(double, double);

  //! Initializes the random number generator and the weights of the net.
  void initWeights(long, double, double);

  //! The core method to set the structure and initialize the bufferes
  void setStructure(Array<int> &);

  //! Based on given connection and weight matrices a network 
  //! is created.
  void setStructure(Array<int> &, Array<double> &);

  //! Replaces the existing weight matrix by a new matrix.
  void setStructure(Array<double> &);

  //! The current structure of the network is replaced by
  //! a new one that is based on the information in file
  //! "filename".
  void setStructure(std::string);

  //! Writes the structure of the network to a file named
  //! "filename".  
  void write(std::string);

  //! Sets the history of the neuron states to certain values. 
  void setHistory(const Array<double> &);

  //! Returns the history of the neuron states.
  Array<double> getHistory();


protected:

  //! The total number of neurons of the network (input, output and hidden). 
  unsigned numberOfNeurons;   

  //! The absolute number of parameters of the network 
  //! (equal to the number of weights)
  unsigned numberOfParameters;


  // the variable that counts the time (backwards!) when processing a data series
  unsigned time;      
  
  // Number of different time steps in the network architecture:
  // delay = 1 means only feed forward, delay > 1 also includes
  // the history of the states into the dynamic.
  unsigned delay;         

  //
  unsigned episode;

  // Number of patterns of the input sequence, which are not
  // considered for evaluating the error and gradient of the
  // error. Nevertheless, these elements are used to modify the
  // internal states of the neurons.
  unsigned warmUpLength;

  // 3-dimensional connection matrix. The first dimension is the time
  // delay, the second and third dimension are the neuron indices of
  // the endpoint and starting point of the connection as given
  // in a "normal" connection matrix.
  ArrayTable<int> connectionMatrix;

  // 3-dimensional weight matrix. The first dimension is the time
  // delay, the second and third dimension are the neuron indices of
  // the endpoint and starting point of the connection, the weight
  // belongs to, as given in a "normal" weight matrix.
  ArrayTable<double> weightMatrix;

  // 1-dimensional array with the number of elements equal to the
  // number of time delays in the structure. The i-th element of this
  // array is set to one, if at least one connection from the i-th memory layer
  // to the feed-forward layer exists. This variable is only used to
  // reduce the computational time.
  ArrayTable<int> delMask;

  // Activation of the neurons prior to the processings, usually equal
  // to the input pattern. Stimulus is a 1-dimensional array with the
  // number of elements equal to the number of neurons.
  ArrayTable<double> stimulus;

  // Activation of the neurons after processing the time series. "Y"
  // is a 2-dimensional array, the first dimension gives the neuron
  // index, the second one the time step counted backwards.
  // Therefore, the second dimension's number of elements is equal to
  // the sum of the length of the time series and the maximum time
  // delay of the structure.
  ArrayTable<double> Y;

  // This array stores the errors of the neurons for every input
  // pattern.
  ArrayTable<double> err;

  // Stores the local delta values which are used to evaluate the
  // gradient.
  ArrayTable<double> delta;

  // Derivative of the error with respect to the weights. This object is 
  // a 3-dimensional array, the first dimension is the time
  // delay, the second and third dimension are the neuron indices of
  // the endpoint and starting point of the connection, the weight
  // belongs to.
  ArrayTable<double> dEw; 

  // Activation function of all neurons.
  virtual double g(double);

  // Computes the derivative of g(a) for all neurons.
  virtual double dg(double);

  // Initializes some internal variables.
  void init0();

  // Writes the values of the weights stored in the weight matrix to 
  // the parameter vector w of ModelInterface.
  void writeParameters();

  // Reads the values of the parameter vector w of ModelInterface
  // and stores these values in the weight matrix.
  void readParameters();

  // Writes the values of the gradient of the error with respect to the
  // different weights from the variable dEw to the variable
  // dedw in the Model Interface.
  void writeGradient();

  // Performs some initializations that are necessary to process the
  // time series.
  void prepareTime(unsigned t);
  
  // Processes a whole time series. After processing the output can be 
  // found in the variable Y.
  void processTimeSeries(const Array<double>& input);

  // Processes one input pattern.
  void processTimeStep();

  // Performs backpropagation through time to calculate the derivative
  // of the error with respect to the weights. The results are stored to
  // dEw.
  void calcGradBPTT();
};

#endif //RNNET_H








