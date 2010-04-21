//===========================================================================
/*!
 *  \file MSERNNet.cpp
 *
 *  \brief Offers the functions to create and to work with a
 *         recurrent neural network.
 *         The network is combined with the mean squared error
 *         measure. This combination is created due to computational
 *         efficiency. 
 *         
 *  \author  M. Toussaint and a ghost
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
 *      $RCSfile: MSERNNet.cpp,v $<BR>
 *      $Revision: 2.3 $<BR>
 *      $Date: 2005/08/03 12:00:15 $
 *
 *  \par Changes:
 *      $Log: MSERNNet.cpp,v $
 *      Revision 2.3  2005/08/03 12:00:15  christian_igel
 *      changes limits/values includes etc.
 *
 *      Revision 2.2  2004/06/01 15:41:44  saviapbe
 *
 *      The bugs in the doxygen's documentation were removed.
 *
 *      Revision 2.1  2004/05/27 15:09:03  saviapbe
 *
 *      Documentation was changed.
 *
 *      Revision 2.0  2003/11/28 16:23:14  shark-admin
 *      Revision tag reset to revision tag 2.x
 *
 *      Revision 1.1.1.1  2003/11/24 13:37:06  shark-admin
 *      INI Administration
 *
 *      Revision 1.2  2002/02/06 15:39:11  rudi
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


#include "ReClaM/MSERNNet.h"

//===========================================================================
/*!
 *  Constructor no. 1
 *
 *  \brief same as #init0();#
 *
 *  \author  not M. Toussaint
 *  \date    2000
 *
 *  \par Changes
 *      none
 *
 *  \par Status
 *      stable
 *
 */  
MSERNNet::MSERNNet(){
  init0();
}


//===========================================================================
/*!
 *  Constructor no. 2
 *
 *  \brief same as #init0(); setStructure(con);#
 *
 *  \param con The 3-dimensional connection matrix used for the
 *             topology of the new network.
 *
 *  \author  not M. Toussaint
 *  \date    2000
 *
 *  \par Changes
 *      none
 *
 *  \par Status
 *      stable
 *
 */  
MSERNNet::MSERNNet(Array<int> con){ 
  init0();
  setStructure(con); 
}


//===========================================================================
/*!
 *  Constructor no. 3
 *
 *  \brief same as #init0(); setStructure(con,wei);#
 *
 *  \param con The 3-dimensional connection matrix used for the
 *             topology of the new network.
 *  \param wei The 3-dimensional weight matrix used for the
 *             topology of the new network.
 *
 *  \author  not M. Toussaint
 *  \date    2000
 *
 *  \par Changes
 *      none
 *
 *  \par Status
 *      stable
 *
 */  
MSERNNet::MSERNNet(Array<int> con,Array<double> wei){
  init0();
  setStructure(con,wei); 
}


//===========================================================================
/*!
 *  Constructor no. 4
 *
 *  \brief same as #init0(); setStructure(filename);#
 *
 *  \param filename Name of the file that contains the information
 *                  for the creation of the network. If the file
 *                  doesn't exist or doesn't have the right format, 
 *                  the method will exit with failure.
 *  \return None.
 *
 *  \author  M. H&uuml;sken
 *  \date    2000
 *
 *  \par Changes
 *      none
 *
 *  \par Status
 *      stable
 *
 */  
MSERNNet::MSERNNet(std::string filename){ 
  init0();
  setStructure(filename); 
}


//===========================================================================
/*!
 *  \brief Sets the length of the warmup sequence.
 *
 *  Usually, when processing a new data series (e.g., by calling
 *  model(...), error(...), or anything like that) then all the
 *  `states' of the network are reset to zero. By `states' I mean the
 *  buffered activations to which time-delayed synapses refer
 *  to. Effectively, this means one assumes a zero activation history.
 *  
 *  The advantage of this is, that it makes the model behavior well
 *  defined. The disadvantage is that you can't predict a time series
 *  well with a zero history. Thus, one should use the first few
 *  inputs of a data series to initialize the network, i.e., to let it
 *  converge into a `normal' dynamic state from which prediction of
 *  new data is possible. This phase is called the warmup
 *  phase. (Internally, the warmup phase differes from subsequent data
 *  only in that it does not contributed to the error functional and
 *  does not induce an error gradient.)
 *  
 *  With this function you can set to time span of the warmup phase
 *  (which also means the amount of data that you `waste' for the
 *  warmup instead of the learning.)
 *
 *  NOTE: Sometimes, e.g., when feeding a sequence of data sequences
 *  to the model, it is desirable not to reset the internal states to
 *  zero. This is the case when you set a WUP to zero.
 *
 *  \param  WUP Length of the warmup sequence, the default value is "0".
 *  \return None.
 *
 *  \author  not M. Toussaint
 *  \date    2000
 *
 *  \par Changes
 *      none
 *
 *  \par Status
 *      stable
 *
 */ 
void MSERNNet::includeWarmUp(unsigned WUP){
  warmUpLength=WUP;
}


//===========================================================================
/*!
 *  \brief Feed a data series to the model.
 *
 *  The output of the network is only stored internally. \em input
 *  must be 2-dimensional, where the first index refers to the time,
 *  the second to the neuron. To view the outputs, please use method
 *  #model(const Array<double>& input,Array<double>& output) instead.
 *  
 *  \param  input Input patterns for the network.
 *  \return None.
 *
 *  \author  M. Toussaint
 *  \date    2000
 *
 *  \par Changes
 *      none
 *
 *  \par Status
 *      stable
 *
 */ 
void MSERNNet::model(const Array<double> &input){
  if(input.ndim()!=2){
    std::cerr << "STOP: MSERNNet::model requires a 2-dim input array\n";
    exit(1);
  }
  if(warmUpLength>0) Y=0.;
  processTimeSeries(input);
}


//===========================================================================
/*!
 *  \brief Feed a data series to the model. The output (i.e., the time
 *  series of activations of the output neurons) is copies into the
 *  #output# buffer.
 *
 *  \param  input  Input patterns for the network.
 *  \param  output Used to store the outputs of the network.
 *  \return None.
 *
 *  \author  M. Toussaint
 *  \date    2000
 *
 *  \par Changes
 *      none
 *
 *  \par Status
 *      stable
 *
 */ 
void MSERNNet::model(const Array<double>& input,Array<double>& output){
  if(output.ndim()!=2){
    std::cerr << "STOP: MSERNNet::model requires a 2-dim output array\n";
    exit(1);
  }
  unsigned t,i,T=output.dim(0),M=output.dim(1);
  model(input);
  for(t=0;t<T;t++) 
    for(i=0;i<M;i++) 
      output(t,i)=Y(numberOfNeurons-1-i,T-1-t);
}



//===========================================================================
/*!
 *  \brief Evaluates the mean squared error of the network output compared 
 *         to given target data.
 *
 *  Given the patterns in \em input and the patterns in \em target
 *  for the comparison to the outputs of the network, the mean
 *  squared error \f$E\f$ as given in the <a href="#_details>details</a>
 *  is calculated.
 *
 *  Keep in mind, that if you have defined a warmup length greater
 *  than zero for the network before, the first \em warmup-length
 *  input patterns are not used for the calculation of the error. Of
 *  course, the number of patterns in \em input must be greater than
 *  the warmup-length otherwise the method will exit with failure.
 *
 *  \param  input  Input patterns for the network.
 *  \param  target Target patterns used for comparison to the outputs
 *                 of the network.
 *  \return The mean squared error \f$E\f$.
 *
 *  \author  M. Toussaint
 *  \date    2000
 *
 *  \par Changes
 *      none
 *
 *  \par Status
 *      stable
 *
 */ 
double MSERNNet::error(const Array<double> &input,const Array<double> &target)
{
  if(target.ndim()!=2){
    std::cerr << "STOP: MSERNNet::error requires a 2-dim target array\n";
    exit(1);
  }

  unsigned i,t,T,M;
  T=target.dim(0);
  M=target.dim(1);

  if(T <= warmUpLength) {
    std::cerr << "number of patterns not larger than warmup length\n";
    exit(1);
  }

  model(input);
  
  err.A=0;
  double totalError=0;
  double z;
  for(t=0;t<warmUpLength;t++)
    for(i=0;i<M;i++){
      err(numberOfNeurons-1-i,T-1-t)=0.0;
    }
  for(t=warmUpLength;t<T;t++)
    for(i=0;i<M;i++){
      z=err(numberOfNeurons-1-i,T-1-t)=Y(numberOfNeurons-1-i,T-1-t)-target(t,i);	
      totalError+=z*z;
    }
  totalError /= (T - warmUpLength) * M; 
  return totalError;
}

//===========================================================================
/*!
 *  \brief Evaluates the error percentage of the network output compared 
 *         to given target data.
 *
 *  Given the patterns in \em input and the patterns in \em target
 *  for the comparison to the outputs of the network, the error
 *  percentage is calculated, i.e. the mean
 *  squared error \f$E\f$ as given in the <a href="#_details>details</a>
 *  is calculated and the error percentage is then given as
 *
 *  \f$
 *      EP = \frac{100}{(outmax - outmin)^2} \ast E
 *  \f$
 *
 *  where \f$outmax\f$ is the maximum value that can be stored
 *  in a \em double value and \f$outmin\f$ is the corresponding
 *  minimum value.
 *  Keep in mind, that if you have defined a warmup length greater than zero 
 *  for the network before, the first \f$warmup-length\f$ input patterns
 *  are not used for the calculation of the error.
 *  Of course, the number of patterns in \em input must be greater
 *  than the warmup-length otherwise the method will exit
 *  with failure.
 *
 *  \param  input  Input patterns for the network.
 *  \param  target Target patterns used for comparison to the outputs
 *                 of the network.
 *  \return The error percentage.
 *
 *  \author  not M. Toussaint
 *  \date    2000
 *
 *  \par Changes
 *      none
 *
 *  \par Status
 *      stable
 *
 */ 
double MSERNNet::errorPercentage(const Array<double> &input,const Array<double> &target)
{
  unsigned t;
#ifndef __SOLARIS_
  double outmax = std::numeric_limits< double >::min( );
  double outmin = std::numeric_limits< double >::max( );
#else
  double outmax = -MAXDOUBLE;
  double outmin = MAXDOUBLE;
#endif

  if(target.ndim() != 2){
	  std::cerr << "STOP: MSERNNet::error requires a 2-dim target array\n";
    exit(EXIT_FAILURE);
  }

  unsigned T, M;
  T = target.dim(0);
  M = target.dim(1);

  if(T <= warmUpLength) {
	  std::cerr << "number of patterns not larger than warmup length\n";
    exit(EXIT_FAILURE);
  }

  model(input);
  
  err.A      = 0;
  double totalError = 0;
  double z;
  for(t = 0; t < warmUpLength; t++)
    for(unsigned i = 0; i < M; i++){
      err(numberOfNeurons-1-i, T-1-t) = 0.0;
    }
  for(t=warmUpLength; t<T; t++)
    for(unsigned i=0; i<M; i++) {
      if(target(t, i) > outmax) outmax = target(t, i);
      if(target(t, i) < outmin) outmin = target(t, i);
      z = err(numberOfNeurons-1-i, T-1-t) = Y(numberOfNeurons-1-i, T-1-t) - target(t, i);	
      totalError += z * z;
    }
  totalError /= (T - warmUpLength) * M; 
  totalError *= 100 / ((outmax - outmin) * (outmax - outmin));
  return totalError;
}


//===========================================================================
/*!
 *  \brief same as #error(input, target);#
 *
 *  \param  input  Input patterns for the network.
 *  \param  target Target patterns used for comparison to the outputs
 *                 of the network.
 *  \return The mean squared error \f$E\f$.
 *
 *  \author  not M. Toussaint
 *  \date    2000
 *
 *  \par Changes
 *      none
 *
 *  \par Status
 *      stable
 *
 */ 
double MSERNNet::meanSquaredError(const Array<double> &input,
				  const Array<double> &target){
  return error(input, target);
}


//===========================================================================
/*!
 *  \brief Evaluates the derivative of the mean squared error.
 *
 *  Given the patterns in \em input and the patterns in \em target
 *  for the comparison to the outputs of the network, the derivatives of
 *  the mean squared error \f$E\f$ as given in the 
 *  <a href="#_details>details</a> are calculated and the results
 *  are stored in ModelInterface::dedw.
 *  Additionally, the mean squared error itself is returned,
 *  if the flag \em returnError is set to "true".
 *  Keep in mind, that if you have defined a warmup length greater than zero 
 *  for the network before, the first \f$warmup-length\f$ input patterns
 *  are not used for the calculation of the error.
 *
 *
 *  \param  input       Input patterns for the network.
 *  \param  target      Target patterns used for comparison to the outputs
 *                      of the network.
 *  \param  returnError Determines whether or not to calculate the error 
 *                      itself. By default the error is calculated.   
 *  \return The error \f$E\f$ if \em returnError is set to "true", "-1"
 *          otherwise.
 *
 *  \author  M. Toussaint
 *  \date    2000
 *
 *  \par Changes
 *      none
 *
 *  \par Status
 *      stable
 *
 */ 
double MSERNNet::derror(const Array<double> &input,
			const Array<double> &target, 
			bool returnError){
  double e;
  e=error(input,target);
  calcGradBPTT();
  writeGradient();
  dedw /= (double) (target.dim(0) - warmUpLength) * target.dim(1);
  if (!returnError)
    e=-1;
  return e;
}  


//===========================================================================
/*!
 *  \brief Returns the connection Matrix.  
 *
 *  The 3-dimensional connection matrix of the network is returned.
 *  The first dimension is used for the number of delays, with
 *  each delay having a 2-dimensional connection matrix as
 *  known from other networks.
 *
 *  \return The 3-dimensional connection matrix.
 *
 *  \author  M. Toussaint
 *  \date    2000
 *
 *  \par Changes
 *      none
 *
 *  \par Status
 *      stable
 *
 */ 
Array<int> & MSERNNet::getConnections(){ 
  return connectionMatrix.A;
}


//===========================================================================
/*!
 *  \brief Returns the weight matrix.
 *
 *  The 3-dimensional weight matrix of the network is returned.
 *  The first dimension is used for the number of delays, with
 *  each delay having a 2-dimensional weight matrix as
 *  known from other networks.
 *
 *  \return The 3-dimensional weight matrix.
 *
 *  \author  M. Toussaint
 *  \date    2000
 *
 *  \par Changes
 *      none
 *
 *  \par Status
 *      stable
 *
 */ 
Array<double> & MSERNNet::getWeights(){ 
  return weightMatrix.A; 
}

//===========================================================================
/*!
 *  \brief Initializes the random number generator and the weights of the net.
 *
 *  same as #Rng::seed(seed); initWeights(min, max);#
 *
 *  \param  seed The initialization value for the random number generator.
 *  \param  min  The minimum possible initialization value.
 *  \param  max  The maximum possible initialization value.
 *  \return None.
 *
 *  \author  not M. Toussaint
 *  \date    2000
 *
 *  \par Changes
 *      none
 *
 *  \par Status
 *      stable
 *
 */ 
void MSERNNet::initWeights(long seed, double min, double max){
  Rng::seed(seed);
  initWeights(min, max);
};


//===========================================================================
/*!
 *  \brief Initializes the weights of the net.
 *
 *  The weights of the network for all delays are initialized 
 *  by uniformally distributed numbers between \em low and \em up.
 *
 *  \param  low The minimum possible initialization value.
 *  \param  up  The maximum possible initialization value.
 *  \return None.
 *
 *  \author  M. Toussaint
 *  \date    2000
 *
 *  \par Changes
 *      none
 *
 *  \par Status
 *      stable
 *
 */ 
void MSERNNet::initWeights(double low,double up){
  unsigned i,j,t;
  for(t=0;t<delay;t++) for(i=0;i<numberOfNeurons;i++) for(j=0;j<numberOfNeurons;j++){
    if(connectionMatrix(t,i,j)) weightMatrix(t,i,j)=Rng::uni(low,up);
    else weightMatrix(t,i,j)=0;
  }
  writeParameters();
}



//===========================================================================
/*!
 *  \brief Based on a given connection matrix a network is created.
 *
 *  This is the core function to set the structure and initialize the
 *  sizes of all fields. All other `setStructure' methods call this
 *  one. Any old values of the network are lost.
 *
 *  One has to pass a 3-dimensional connection matrix \em mat to this
 *  function. The first dimension indicates the time delay of the
 *  respective connection. E.g., mat[0] is the ordinary 2-dimensional
 *  feedforward (triangular) connection matrix (where the connections
 *  have no time delay), mat[1] is the (potentially fully connected)
 *  connection matrix for connections of time delay 1, etc.
 *
 *  It is perfectly ok, if some of these `layers' mat[t] are
 *  completely zero: the function will detect this (and store it as a
 *  0 in the delMask(t)) and not waste time for this layer when
 *  processing the network. 
 *
 *  Notice, that the weight matrix #w# of the network is only adapted
 *  to the new size, but the matrix is not explicitly initialized with
 *  zero weights.
 *
 *  \param  mat The 3-dimensional connection matrix.

 *  \return None.
 *
 *  \author  M. Toussaint
 *  \date    2000
 *
 *  \par Changes
 *      none
 *
 *  \par Status
 *      stable
 *
 */ 
void MSERNNet::setStructure(Array<int> &mat){
  unsigned i,j,t;
  
  delay=mat.dim(0);
  numberOfNeurons=mat.dim(1);
  weightMatrix.resize(delay,numberOfNeurons,numberOfNeurons,false);
  dEw.resize(delay,numberOfNeurons,numberOfNeurons,false);
  stimulus.resize(numberOfNeurons,false);
  connectionMatrix.resize(delay,numberOfNeurons,numberOfNeurons,false);
  delMask.resize(delay,false);
  
  // ergaenzt 5.9. huesken
  // M.Toussaint: wieso denn das??? Y. wird erst in prepareTime initialisiert...
  Y.resize(numberOfNeurons,delay,false);
  Y.A = 0;
  delta.resize(numberOfNeurons,delay,false);
  delta.A = 0;
  
  connectionMatrix.A=mat;
  numberOfParameters=0;
  for(t=0;t<delay;t++){
    delMask(t)=0;
    for(i=0;i<numberOfNeurons;i++) for(j=0;j<numberOfNeurons;j++) if(connectionMatrix(t,i,j)){
      numberOfParameters++;
      delMask(t)=1;
    }
  }
  w.resize(numberOfParameters,false);
  dedw.resize(numberOfParameters,false);
}



//===========================================================================
/*!
 *  \brief Based on given connection and weight matrices a network 
 *         is created.
 *
 *  same as #setStructure(con); setStructure(wei);#
 *
 *  \param  con The 3-dimensional connection matrix, where the
 *              first dimension is used for the different time
 *              steps with each time step having a 2-dimensional
 *              #connectionMatrix.
 *  \param  wei The new 3-dimensional weight matrix, where the
 *              first dimension is used for the different time
 *              steps with each time step having a 2-dimensional
 *              #weightMatrix.
 *  \return None.
 *
 *  \author  not M. Toussaint
 *  \date    2000
 *
 *  \par Changes
 *      none
 *
 *  \par Status
 *      stable
 *
 */ 
void MSERNNet::setStructure(Array<int> &con, Array<double> &wei){
  setStructure(con);
  setStructure(wei);
}


//===========================================================================
/*!
 *  \brief Replaces the existing weight matrix by a new matrix.
 *
 *  The current weight matrix is replaced by the new weight
 *  matrix \em wei. 
 *  This means NOT, that the whole structure of the network
 *  is changed. The new weight matrix must be compatible
 *  to the existing connection matrix.
 *
 *  \param  wei The new 3-dimensional weight matrix, where the
 *              first dimension is used for the different time
 *              steps with each time step having a 2-dimensional
 *              #weightMatrix.
 *  \return None.
 *
 *  \author  not M. Toussaint
 *  \date    2000
 *
 *  \par Changes
 *      none
 *
 *  \par Status
 *      stable
 *
 */ 
void MSERNNet::setStructure(Array<double> &wei){ 
  weightMatrix.A=wei; 
  writeParameters();
}



//===========================================================================
/*!
 *  \brief The current structure of the network is replaced by
 *         a new one that is based on the information in file
 *         "filename".
 *
 *  A file is used to replace the current structure of the network. 
 *  This file must have the following content:
 *  The first line contains the number of time steps used,
 *  followed by \f$number\ of\ time\ steps\f$ connection matrices,
 *  followed by \f$number\ of\ time\ steps\f$ weight matrices.
 *
 *  \param filename Name of the file that contains the information
 *                  for the modification of the structure.
 *  \return None.
 *
 *  \author  M. H&uuml;sken
 *  \date    2000
 *
 *  \par Changes
 *      none
 *
 *  \par Status
 *      stable
 *
 */  
void MSERNNet::setStructure(std::string filename){
  // lade Netzdaten
  std::ifstream is;
  is.open(filename.c_str());
  unsigned nDelay;
  is >> nDelay;
  
  std::vector<int> l;
  const unsigned Len = 256;
  char buffer[Len];
  unsigned n = 0;
  is >> std::setw(Len) >> buffer;

  while(is) {
    n++;
    l.push_back(atoi(buffer));
    if(is.peek()=='\n') {
      break;
    }

    is >> std::setw(Len) >> buffer;
  }
  
  if(n == 0) {
	std::cerr << "no matrices given in net file";
    exit(EXIT_FAILURE);
  }
  
  Array<int> con (nDelay,n,n);    
  Array<double> wei (nDelay,n,n);
  unsigned i, j, k;
  
  for (i=0; i<con.dim(0); i++)
    for (j=0; j<con.dim(1); j++)
      if ((i==0) && (j==0)){
	for(k = 0; k < n; k++) con(i, j, k) = l[k];
      }
      else
	for (k=0; k<con.dim(2); k++){
	  is >> con(i, j, k);
	  if(!is) {
		std::cerr << "not enough entries to build connection matrix" << std::endl;
	    exit(EXIT_FAILURE);
	  }
	}
  
  for (i=0; i<wei.dim(0); i++)
    for (j=0; j<wei.dim(1); j++)
      for (k=0; k<wei.dim(2); k++){
	is >> wei(i, j, k);
	if(!is) {
	  std::cerr << "not enough entries to build weight matrix" << std::endl;
	  exit(EXIT_FAILURE);
	}
      }
  
  setStructure(con,wei);
}


//===========================================================================
/*!
 *  \brief Writes the structure of the network to a file named
 *         "filename".  
 *
 *  The structure of the network is written in the following
 *  format:
 *
 *  The first line only contains a single value, that denotes
 *  the number of time steps used.
 *  Following are \f$number\ of\ time\ steps\f$ connection
 *  matrices, followed by \f$number\ of\ time\ steps\f$
 *  weight matrices.
 *
 *  \param  filename Name of the file, the network structure is
 *                   written to.
 *  \return None.
 *
 *  \author  not M. Toussaint
 *  \date    2000
 *
 *  \par Changes
 *      none
 *
 *  \par Status
 *      stable
 *
 */ 
void MSERNNet::write(std::string filename){
  unsigned i,j,k;
  FILE *fp=fopen(filename.c_str(),"wt");
  
  Array<int> con = getConnections();
  unsigned nDelay = con.dim(0);
  fprintf(fp,"%u\n",nDelay);
  
  for (k=0; k<nDelay; k++){  
    for (i=0; i<con.dim(1); i++){
      for (j=0; j<con.dim(2); j++){
	fprintf(fp,"%1d",con(k,i,j));
	if (j<con.dim(2)-1) fprintf(fp," ");
      }
      fprintf(fp,"\n");
    }
    fprintf(fp,"\n");
  }
  
  Array<double> wei = getWeights();
  for (k=0; k<nDelay; k++){  
    for (i=0; i<wei.dim(1); i++){
      for (j=0; j<wei.dim(2); j++){
	fprintf(fp,"%10.8e",wei(k,i,j));
	if (j<wei.dim(2)-1) fprintf(fp," ");
      }
      fprintf(fp,"\n");
    }
    fprintf(fp,"\n");
  }
  fclose(fp);
}


//===========================================================================
/*!
 *  \brief Sets the history of the neuron states to certain values.
 *
 *  In case you chose no warmup phase (warmup length set to zero)
 *  then, the network will not be reinitialized when new data is
 *  feeded in. This might make sense if previous data was processes
 *  and the current state of the model is well-defined and can ne used
 *  for further predictions.
 *
 *  Alternatively, with this function, you can explicitly set the
 *  history of the network -- I guess this makes only sense when you
 *  have recorded a history earlier or if you set it explicitly to
 *  zero.
 *
 *  The first index of the array corresponds to the index of the
 *  neuron (and must have dimension equal to the total number of
 *  neurons). The second index corresponds to the time IN REVERSE
 *  ORDER (the index denotes `the time before now') and it must have
 *  the dimensionality #delay#.
 *
 *  \param  Ystate The new values for the history.
 *  \return None.
 *
 *  \author  M. Toussaint
 *  \date    2000
 *
 *  \par Changes
 *      none
 *
 *  \par Status
 *      stable
 *
 */ 
void MSERNNet::setHistory(const Array<double> &Ystate)
{
  if(Ystate.dim(0)!=numberOfNeurons || Ystate.dim(1)!=delay){
	std::cerr << "STOP: MSERNNet::setHistory: wrong dimensions\n";
    exit(1);
  }
  unsigned i,j;
  for(i=0;i<numberOfNeurons;i++) 
    for(j=0;j<delay;j++) 
      Y(i,j)=Ystate(i,j);
}



//===========================================================================
/*!
 *  \brief Returns the history of the neuron states.
 *
 *  This returns the current history of neuron activations. Use it for
 *  #setHistory#. See the explanations of the history format there.
 *
 *  \return The current content of the history.
 *
 *  \author  M. Toussaint
 *  \date    2000
 *
 *  \par Changes
 *      none
 *
 *  \par Status
 *      stable
 *
 */ 
Array<double> MSERNNet::getHistory(){
  unsigned i,j;
  Array<double> Ystate(numberOfNeurons,delay);
  for(i=0;i<numberOfNeurons;i++)
    for(j=0;j<delay;j++)
      Ystate(i,j)=Y(i,j);
  return Ystate;
}


//###########################################################################
//
// Private Methods:
//
//###########################################################################


// Activation function of all neurons.
double MSERNNet::g(double a){ 
  return 1/(1+exp(-a));
}


// Computes the derivative of g(a) for all neurons.
double MSERNNet::dg(double ga){
  return ga*(1-ga); 
}


// Initializes some internal variables.
void MSERNNet::init0(){
  time=numberOfNeurons=delay=episode=numberOfParameters=0;
  warmUpLength=0;
  Y.resize(0,0,false);
  err.resize(0,0,false);
  delta.resize(0,0,false);
}


// Writes the values of the weights stored in the weight matrix to 
// the parameter vector w of ModelInterface.
void MSERNNet::writeParameters(){
  unsigned i,j,t,n=0;
  for(t=0;t<delay;t++) for(i=0;i<numberOfNeurons;i++) for(j=0;j<numberOfNeurons;j++)
    if(connectionMatrix(t,i,j))	w(n++)=weightMatrix(t,i,j);
}


// Reads the values of the parameter vector w of ModelInterface
// and stores these values in the weight matrix.
void MSERNNet::readParameters(){
  unsigned i,j,t,n=0;
  for(t=0;t<delay;t++) for(i=0;i<numberOfNeurons;i++) for(j=0;j<numberOfNeurons;j++)
    if(connectionMatrix(t,i,j)) weightMatrix(t,i,j)=w(n++); else weightMatrix(t,i,j)=0;
}


// Writes the values of the gradient of the error with respect to the
// different weights from the variable dEw to the variable
// dedw in the Model Interface. 
void MSERNNet::writeGradient(){
  unsigned i,j,t,n=0;
  for(t=0;t<delay;t++) for(i=0;i<numberOfNeurons;i++) for(j=0;j<numberOfNeurons;j++)
    if(connectionMatrix(t,i,j)) dedw(n++)=dEw(t,i,j);
}


// Performs some initializations that are necessary to process the
// time series.
void MSERNNet::prepareTime(unsigned t){
  unsigned i,j;
  ArrayTable<double> Ystate,Dstate;
  Ystate.resize(numberOfNeurons,delay,false);
  Dstate.resize(numberOfNeurons,delay,false);
  for(i=0;i<numberOfNeurons;i++) for(j=0;j<delay;j++){
    Ystate(i,j)=Y(i,j);
    Dstate(i,j)=delta(i,j);
  }
  if(t!=episode){
    episode=t;
    Y.resize(numberOfNeurons,episode+delay,false);
    Y.A = 0;
    delta.resize(numberOfNeurons,episode+delay,false); 
    delta.A = 0;
    err.resize(numberOfNeurons,episode+delay,false);
    err.A = 0;
  }
  for(i=0;i<numberOfNeurons;i++) for(j=0;j<delay;j++){
    Y(i,t+j)=Ystate(i,j);
    delta(i,t+j)=Dstate(i,j);
  }
  time=t;
}


// Processes a whole time series. After processing the output can be 
// found in the variable Y.
void MSERNNet::processTimeSeries(const Array<double>& input){
  unsigned t,i;
  if(input.dim(1)>=numberOfNeurons){
	std::cerr << "STOP: structure has not enough neurons to handle input"
	 << " - use setStructure(...)\n";
    exit(1);
  }
  
  readParameters();
  prepareTime(input.dim(0));
  stimulus.A=0;
  for(t=0;t<input.dim(0);t++){
    for(i=0;i<input.dim(1);i++) stimulus(i)=input(t,i);
    processTimeStep();
  }
}


// Processes one input pattern.
void MSERNNet::processTimeStep(){
  unsigned i,j,t;
  double z;
  if(!time){ 
    fprintf(stderr,"MSERNNet::process: no time prepared!"); 
    exit(0); 
  }
  time--;
  for(i=0;i<numberOfNeurons;i++){
    z=stimulus(i);
    z+=weightMatrix(0,i,i);
    if(delMask(0)) 
      for(j=0;j<i;j++)
	z+=weightMatrix(0,i,j)*Y(j,time);
    for(t=1;t<delay;t++)
      if(delMask(t))
	for(j=0;j<numberOfNeurons;j++)
	  z+=weightMatrix(t,i,j)*Y(j,time+t);
    Y(i,time)=g(z);
  }
}


// Performs backpropagation through time to calculate the derivative
// of the error with respect to the weights. The results are stored to
// dEw.
void MSERNNet::calcGradBPTT(){
  int i,j,t,s,n;
  double z;
  
  for(t=0;t<(int)episode;t++) for(i=numberOfNeurons-1;i>=0;i--){
    z=err(i,t);
    if(delMask(0)) for(j=numberOfNeurons-1;j>i;j--)
      z+=delta(j,t)*weightMatrix(0,j,i);
    for(s=1;s<(int)delay && s<=t;s++) if(delMask(s)) for(j=numberOfNeurons-1;j>=0;j--)
      z+=delta(j,t-s)*weightMatrix(s,j,i);
    delta(i,t)=dg(Y(i,t)) * z;
  }
  
  for(t=0;t<(int)delay;t++) if(delMask(t)) for(i=0;i<(int)numberOfNeurons;i++){
    if (t>0) n=numberOfNeurons; else n=i;
    for(j=0;j<n;j++) if(connectionMatrix(t,i,j)){
      dEw(t,i,j)=0;
      for(s=0;s+t<(int)episode;s++) dEw(t,i,j)+=delta(i,s) * Y(j,t+s);
      dEw(t,i,j)*=2.0;
    }
  }
  
  for(i=0;i<(int)numberOfNeurons;i++) if(connectionMatrix(0,i,i)){
    dEw(0,i,i)=0;
    for(s=0;s<(int)episode;s++) dEw(0,i,i)+=delta(i,s);
    dEw(0,i,i)*=2.0;
  }
}
