//===========================================================================
/*!
 *  \file DF_CrossEntropyIndependent.h
 *
 *  \brief Error measure for classification tasks of non exclusive attributes
 *         that can be used for model training. 
 * 
 *  If your model should return a vector whose components are connected to 
 *  multiple mutually independent attributes and the \em k-th component of the 
 *  output vector is representing the probability of the presence of the \em k-th 
 *  attribute given any input vector, 'DF_CrossEntropyIndependent' is the adequate 
 *  error measure for model-training. For \em C>1, dimension of model's output 
 *  and every output dimension represents a single attribute or class respectively,
 *  it follows the formular
 *  \f[
 *      E = - \sum_{i=1}^N \sum_{k=1}^{C} \{tar^i_k \ln model_k(in^i) + (1-tar^i_k) \ln 
 *          (1-model_k(in^i))\}
 *  \f] 
 *  where \em i runs over all input patterns.
 *  This error functional can be derivated and so used for training. In case 
 *  of only one single output dimension 'DF_CrossEntropyIndependent' returns actually the
 *  true cross entropy for two classes, using the formalism
 *  \f[
 *      E = - \sum_{i=1}^N \{tar^i \ln model(in^i) + (1-tar^i) \ln 
 *          (1-model(in^i))\}
 *  \f]
 *  For theoretical reasons it is suggested to use for neural networks
 *  the logistic sigmoid activation function at the output neurons. 
 *  However, actually every sigmoid activation could be applied to the output 
 *  neurons as far as the image of this function is identical to the intervall
 *  \em [0,1]. 
 *  In this implementation every target value to be chosen from {0,1} 
 *  (binary encoding). For detailed information refer to 
 *  (C.M. Bishop, Neural Networks for Pattern Recognition, Clarendon Press 1996, Chapter 6.8.)        
 *
 *  This implementation of 'CrossEntropyIndependent' performs more efficient than the implemetation given
 *  in 'CrossEntropyIndependent.h' for more than one output dimensions, because redundant calculations of 
 *  the outer derivatives are circumvented.
 * 
 *  \author  S. Wiegand
 *  \date    2003
 *
 *  \par Copyright (c) 1999-2001:
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
 *  \par File and Revision:
 *      $RCSfile: DF_CrossEntropyIndependent.h,v $<BR>
 *      $Revision: 2.1 $<BR>
 *      $Date: 2004/06/01 15:41:45 $
 *
 *  \par Changes:
 *      $Log: DF_CrossEntropyIndependent.h,v $
 *      Revision 2.1  2004/06/01 15:41:45  saviapbe
 *
 *      The bugs in the doxygen's documentation were removed.
 *
 *      Revision 2.0  2003/11/28 16:23:22  shark-admin
 *      Revision tag reset to revision tag 2.x
 *
 *      Revision 1.1.1.1  2003/11/24 13:37:06  shark-admin
 *      INI Administration
 *
 *
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

 */
//===========================================================================


#ifndef DF_CROSS_ENTROPY_INDEPENDENT_H
#define DF_CROSS_ENTROPY_INDEPENDENT_H

#include <cmath>
#include "ReClaM/ModelInterface.h"

//===========================================================================
/*!
  *  \brief Error measure for classification tasks of non exclusive attributes
 *         that can be used for model training. 
 * 
 *  If your model should return a vector whose components are connected to 
 *  multiple mutually independent attributes and the \em k-th component of the 
 *  output vector is representing the probability of the presence of the \em k-th 
 *  attribute given any input vector, 'DF_CrossEntropyIndependent' is the adequate 
 *  error measure for model-training. For \em C>1, dimension of model's output 
 *  and every output dimension represents a single attribute or class respectively,
 *  it follows the formular
 *  \f[
 *      E = - \sum_{i=1}^N \sum_{k=1}^{C} \{tar^i_k \ln model_k(in^i) + (1-tar^i_k) \ln 
 *          (1-model_k(in^i))\}
 *  \f] 
 *  where \em i runs over all input patterns.
 *  This error functional can be derivated and so used for training. In case 
 *  of only one single output dimension 'DF_CrossEntropyIndependent' returns actually the
 *  true cross entropy for two classes, using the formalism
 *  \f[
 *      E = - \sum_{i=1}^N \{tar^i \ln model(in^i) + (1-tar^i) \ln 
 *          (1-model(in^i))\}
 *  \f]
 *  For theoretical reasons it is suggested to use for neural networks
 *  the logistic sigmoid activation function at the output neurons. 
 *  However, actually every sigmoid activation could be applied to the output 
 *  neurons as far as the image of this function is identical to the intervall
 *  \em [0,1]. 
 *  In this implementation every target value to be chosen from {0,1} 
 *  (binary encoding). For detailed information refer to 
 *  (C.M. Bishop, Neural Networks for Pattern Recognition, Clarendon Press 1996, Chapter 6.8.)        
 *
 *  This implementation of 'CrossEntropyIndependent' performs more efficient than the implemetation given
 *  in 'CrossEntropyIndependent.h' for more than one output dimensions, because redundant calculations of 
 *  the outer derivatives are circumvented.        
 *
 *
 *  \author  S. Wiegand
 *  \date    2003
 *
 *  \par Changes:
 *       none
 *
 *  \par Status:
 *      stable
 */
class DF_CrossEntropyIndependent : virtual public ModelInterface {
    // This value is used instead of log of a small number to avoid numerical instabilities. 
  // This value is chosen, as minLog = log(0) on a unix system.
  const static double minLog=-730;
  const static double maxDeriv=1e300;
 public:

//===========================================================================
/*!
 *  \brief Calculates the CrossEntropyIndependent error.
 *
 *  Calculates the CrossEntropyIndependent error function for \em N patterns and 
 *  \em C>1 independent attributes within the output vector via
 *  \f[
 *      E = - \sum_{i=1}^N \sum_{k=1}^{C} \{tar^i_k \ln model_k(in^i) + (1-tar^i_k) \ln 
 *          (1-model_k(in^i))\}
 *  \f] 
 *  or the true cross entropy for only one single output dimension via
 *  \f[
 *      E = - \sum_{i=1}^N \{tar^i \ln model(in^i) + (1-tar^i) \ln 
 *          (1-model(in^i))\}
 *  \f]
 *      \param  in Input vector for the model.
 *      \param  target Target vector.
 *      \return The error \em E.
 *
 *  \author  S. Wiegand
 *  \date    2003
 *
 *  \par Changes
 *       none
 *
 *  \par Status
 *      stable
 *
 */  
  double error(const Array<double> &in, const Array<double> &target) {
    double ce=0;
    if (outputDimension>1)
      // more then one output neuron
      if(in.ndim() == 1) {
	Array<double> output(target.dim(0));
	model(in, output);
	for(unsigned c = 0; c < target.dim(0); c++)
	  if((output(c)!=0) && (1.0-output(c)!=0))
	    ce -= target(c) * log( output(c) ) + ( 1. - target(c) ) * log( 1. - output(c) );
	  else if(output(c)!=target(c))
	    ce -= minLog;
      } else {
	Array<double> output(target.dim(1));
	for(unsigned pattern = 0; pattern < in.dim(0); ++pattern) {
	  model(in[pattern], output);
	  //	  writeArray(output,std::cout);
	  for(unsigned c = 0; c < target.dim(1); c++)
	    if((output(c)!=0) && (1.0-output(c)!=0))
		ce -= target(pattern, c) * log( output(c) ) + (1. - target(pattern, c)) * log( 1. - output(c) );
	    else if(output(c)!=target(pattern, c))
		ce -= minLog;	
	}
      }
    else{
      // only one output neuron
      if(in.ndim() == 1) {
	// one pattern
	Array<double> output(1);
	model(in, output);
	if ((output(0)!=0) && (1.0-output(0)!=0))
	  ce -= (target(0)*log(output(0))+(1.0-target(0))*log(1.0-output(0)));
	else if (output(0)!=target(0))
	  ce -= minLog;
      } else {
	// more then one pattern
	Array<double> output(1);
	for(unsigned pattern = 0; pattern < in.dim(0); ++pattern) {
	  model(in[pattern], output);
	  if ((output(0)!=0) && (1.0-output(0)!=0))
	    ce -= target(pattern,0)*log(output(0))+(1.0-target(pattern,0))*log(1.0-output(0));
	  else if (output(0)!=target(pattern,0))
	    ce -= minLog;	  
	}
      }
    }
    return ce;
  }
  
//===========================================================================
/*!
 *  \brief Calculates the derivatives of the CrossEntropyIndependent error 
 *         (see #error) with respect to the parameters ModelInterface::w.
 *  
 *  Calculates the CrossEntropyIndependent error derivative for \em N patterns and 
 *  \em C>1 independent attributes within the output vector corresponding to model parameter 
 *  \em w via
 *  \f[
 *      \frac{\partial E}{\partial w} = - \sum_{i=1}^N \sum_{k=1}^{C} 
 *                                        \left\{ \frac{tar^i_k - model_k(in^i)}{model_k(in^i)\cdot(1-model_k(in^i))}
 *					  \cdot\frac{\partial model_k(in^i)}{\partial w}\right\}
 *  \f] 
 *  or the derivative of the true cross entropy for only one single output dimension via
 *  \f[
 *  \frac{\partial E}{\partial w} = - \sum_{i=1}^N 
 *                                    \left\{ \frac{tar^i - model(in^i)}{model(in^i)\cdot(1-model(in^i))}
 *				      \cdot\frac{\partial model(in^i)}{\partial w}\right\}
 *  \f]
 *  Usually, as a byproduct of the calculation of the derivative one gets the 
 *  CrossEntropyIndependent error itself very efficiently. Therefore, the method #error
 *  gives back this value. This additional effect can be switched of by means
 *  of the third parameter (returnError = false).
 *
 *      \param  in          Input vector for the model.
 *      \param  target      Target vector.
 *      \param  returnError Determines whether or not to calculate the cross
 *                          entropy error itself. By default the
 *                          error is calculated.  
 *      \return The CrossEntropyIndependent error if \em returnError is set to "true", 
 *              "-1" otherwise.
 *
 *  \author  S. Wiegand
 *  \date    2003
 *
 *  \par Changes
 *       none
 *
 *  \par Status
 *      stable
 *
 */  
  double derror(const Array<double> &in, const Array<double> &target, bool returnError=true) {
    dedw = 0;
    double ce = 0;
    coeffs.resize(outputDimension, false);
    coeffs = 0;
    if (outputDimension>1)
      // more then one output neuron
      if(in.ndim() == 1) {
	// one input pattern
 	Array<double> output(target.dim(0));
	model(in, output);
	if (returnError){
	  // calculate CrossEntropyIndependent
	  for(unsigned c = 0; c < target.dim(0); c++)
	    if((output(c)!=0) && (1.0-output(c)!=0))
	      ce -= target(c) * log( output(c) ) + ( 1. - target(c) ) * log( 1. - output(c) );
	    else if(output(c)!=target(c))
	      ce -= minLog;
	}
	// calculate derivative dedw
	  for(unsigned c = 0; c < output.nelem(); c++)
	    if((output(c)!=0) && (1.0-output(c)!=0))
	      coeffs(c) = (  (target(c) - output(c))/( (1. - output(c)) * output(c) )  );
	    else if(output(c)!=target(c)) 
	      coeffs(c) = ( 2.0*target(c)-1.0 ) * maxDeriv;
	  df(coeffs, dedw);
      } else {
	// more than one input pattern
	Array<double> output(target.dim(1));
	for(unsigned pattern = 0; pattern < in.dim(0); ++pattern) {
	  model(in[pattern], output);
	  if (returnError){
	    // calculate cross-entropy
	    for(unsigned c = 0; c < target.dim(1); c++)
	      if((output(c)!=0) && (1.0-output(c)!=0))
		ce -= target(pattern, c) * log( output(c) ) + (1. - target(pattern, c)) * log( 1. - output(c) );
	      else if(output(c)!=target(pattern, c))
		ce -= minLog;	
	  }
	  // calculate dericative dedw
	  for(unsigned c = 0; c < output.nelem(); c++)
	    if((output(c)!=0) && (1.0-output(c)!=0))
	      coeffs(c) =  (  (target(pattern,c) - output(c))/( (1. - output(c)) * output(c) )  );
	    else if(output(c)!=target(pattern, c))
	      // the sign changes due to the value of the target, maxDeriv is devided by reason of computational bounds
	      coeffs(c) = ( 2.0*target(pattern,c)-1.0 ) * maxDeriv/in.dim(0);
	  df(coeffs,dedw); 
	}
      }
    else {
      // one output neuron
      if(in.ndim() == 1) {
	// only one pattern
	Array<double> output(1);
	model(in, output);
	if (returnError){
	  // calculate ce
	  if ((output(0)!=0) && (1. - output(0)!=0))
	    ce -= (target(0) * log(output(0)) + (1. -target(0)) * log(1. - output(0)));
	  else if (output(0) != target(0)) 
	    ce -= minLog;
	}
	// calculate derivative dedw
	if ((output(0)!=0) && (1. - output(0) != 0)){
	  coeffs(0) = (target(0) - output(0)) / (output(0) * (1. - output(0)));
	  df(coeffs, dedw);
	} else { 
	  if (output(0) != target(0)){
	    // the sign changes due to the value of the target,
	    coeffs(0) = (2. * target(0) - 1.) * maxDeriv;
	    df(coeffs, dedw);
	  }
	}
      } else {
	// a number of patterns
	Array<double> output(1);
	for(unsigned pattern = 0; pattern < in.dim(0); ++pattern) {
	  model(in[pattern], output);//
	  if (returnError){
	    // calculate ce
	    if ((output(0) != 0) && (1. - output(0) != 0))
	      ce -= target(pattern, 0) * log(output(0)) + (1. - target(pattern, 0)) * log(1. - output(0));
	    else  if (output(0)!=target(pattern,0))
	      ce -= minLog;
	  }
	  // calculate derivative dedw
	  if ((output(0)!=0) && (1. - output(0)!=0)){
	    coeffs(0)=(target(pattern, 0)-output(0)) / (output(0) * (1. - output(0)));
	    df(coeffs, dedw);
	  } else { 
	    if (output(0) != target(pattern, 0)){
	      // the sign changes due to the value of the target, maxDeriv is devided by reason of computational bounds
	      coeffs(0)=(2. * target(pattern, 0) - 1.) * maxDeriv / in.dim(0);
	      df(coeffs, dedw);
	    }
	  }
	}
      }
    }
    if (!returnError)
      ce = -1;
    return ce;
  }
  
  private:

  Array<double> coeffs;  
};
#endif








