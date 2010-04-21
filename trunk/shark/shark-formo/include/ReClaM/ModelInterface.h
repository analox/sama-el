
#ifndef _NEW_RECLAM_VERSION_


//===========================================================================
/*!
 *  \file ModelInterface.h
 *
 *  \brief Realizes the communication between the different modules.
 *
 *  To provide flexibility ReClaM offers different modules that can
 *  be put together to form an environment for solving regression
 *  and classification tasks.<BR>
 *  The choice of three module types is necessary: A parametric model,
 *  an error function module and an optimization algorithm model.
 *  All modules exchange parameters and other information, so the task
 *  of this class is to offer this communication.
 *
 *  \author  C. Igel
 *  \date    1999
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
 *      $RCSfile: ModelInterface.h,v $<BR>
 *      $Revision: 2.2 $<BR>
 *      $Date: 2005/12/30 10:50:17 $
 *
 *  \par Changes:
 *      $Log: ModelInterface.h,v $
 *      Revision 2.2  2005/12/30 10:50:17  glasmtbl
 *      *** empty log message ***
 *
 *      Revision 2.0  2003/11/28 16:23:22  shark-admin
 *      Revision tag reset to revision tag 2.x
 *
 *      Revision 1.1.1.1  2003/11/24 13:37:06  shark-admin
 *      INI Administration
 *
 *      Revision 1.6  2002/09/24 08:27:03  rudi
 *      doxygen comments adapted to the new version; new method "df" added.
 *
 *      Revision 1.5  2002/05/16 13:26:03  rudi
 *      doxygen commands added/modified.
 *
 *      Revision 1.4  2001/11/30 14:58:09  rudi
 *      de -> dedw
 *      dw -> dmdw
 *      w, dedw, dmdw protected
 *      fderror removed
 *      derror returns error, third parameter
 *      fdmodel(.,.) -> dmodel(.,.)
 *      regularization removed
 *      doxygen comments added.
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
 *
 */
//===========================================================================

#ifndef MODEL_INTERFACE_H 
#define MODEL_INTERFACE_H 

#include <iostream>
#include "Array/ArrayIo.h"

//===========================================================================
/*!
 *  \brief Realizes the communication between the different modules.
 *
 *  To provide flexibility ReClaM offers different modules that can
 *  be put together to form an environment for solving regression
 *  and classification tasks.<BR>
 *  The choice of three module types is necessary: A parametric model,
 *  an error function module and an optimization algorithm model.
 *  All modules exchange parameters and other information, so the task
 *  of this class is to offer this communication.<BR>
 *  All of the methods of the class are virtual and should to be 
 *  implemented by the different modules.
 *
 *  \author  C. Igel
 *  \date    1999
 *
 *  \par Changes:
 *      none
 *
 *  \par Status:
 *      stable
 */
class ModelInterface{
	public:
  
//===========================================================================
/*!
	 *  \brief Constructs an object with a default value for #epsilon.
				*
				*  This constructor is defined in order to prevent the implicit definition
				*  by the compiler.
				*
				*  \author  C. Igel
				*  \date    1999
				*
				*  \par Changes
				*      none
				*
				*  \par Status
				*      stable
				*
 */
					   ModelInterface(double e = 1e-8) : epsilon(e) 
			   {
			   };

//===========================================================================
/*!
				*  \brief Destructs an object.
						   *
						   *  Don't cause exceptions in destructors.
						   *
						   *  Because destructors are called in the process of throwing other
						   *  exceptions, you'll never want to throw an exception in a destructor
						   *  or cause another exception to be thrown by some action you perform in
						   *  the destructor. If this happens, it means that a new exception may be
						   *  thrown before the catch-clause for an existing exception is reached,
						   *  which will cause a call to terminate().
						   *
						   *  This means that if you call any functions inside a destructor that
						   *  may throw exceptions, those calls should be within a try block in the
						   *  destructor, and the destructor must handle all exceptions itself.
						   *  None must escape from the destructor.
						   *
						   *  \author  C. Igel
						   *  \date    1999
						   *
						   *  \par Changes
						   *      none
						   *
						   *  \par Status
						   *      stable
						   *
 */
								  virtual ~ModelInterface() 
						  {
						  };

//===========================================================================
/*!
						   *  \brief Sets the value of #epsilon.
									  *
									  *  The value of #epsilon is utilized for a numeric estimation of the
									  *  derivative of the model and error, respectively, with respect to the 
									  *  weights. #epsilon should be choosen as a small value. The default value
									  *  is #epsilon = 1e-8.
									  *
									  *      \param  e New value for epsilon.
									  *      \return None.
									  *
									  *  \author  C. Igel
									  *  \date    1999
									  *
									  *  \par Changes
									  *      none
									  *
									  *  \par Status
									  *      stable
									  *
 */
											 void setEpsilon (double e) 
									 { 
										 epsilon = e; 
									 };


									 //
  // Interface variables
									 //
 
	protected:
  //! Parameter vector with weight and bias values.
  //! #w is a 1-dimensional array with the number of elements equal 
  //! to the number of free parameters of the model (i.e., equal 
  //! to the number of weights of a neural network).
		Array<double> w;

  //! Derivatives of the model outputs with respect to the parameters.
  //! #dmdw is a two-dimensional array, dmdw(i,j) equals the derivative
  //! of te i-th output neuron with respect to the j-th weight.
		Array<double> dmdw; 

  //! Derivatives of the error with respect to the parameters. This
  //! array is one-dimensional, the i-th element contains the derivative
  //! of the error with respect to the i-th weight.
		Array<double> dedw;

	public:
  //! Dimension of input vector x.
		unsigned      inputDimension;

  //! Dimension of output vector y.
		unsigned      outputDimension;

//===========================================================================
/*!
		 *  \brief Error of the model. 
					*
					*  This method calculates the error of the model. This function is 
					*  usually specified by classes that define particular error measures.
					*  It is one of the basic needs to specify this method to make ReClaM run.
					*  #target is compared with the model's output for stimulus #input.
					*
					*      \param  input Vector of input values.
					*      \param  target Vector of output values.
					*      \return The error of the model.
					*
					*  \author  C. Igel
					*  \date    1999
					*
					*  \par Changes
					*      none
					*
					*  \par Status
					*      stable
					*
 */ 

						   virtual double error   (const Array<double> &input, const Array<double> &target) = 0;
  

//===========================================================================
/*!
					*  \brief Calculates the field #dedw, the derivative of the error
							   *  \f$E\f$ with respect to the parameters #w.
							   *
							   *  Calculates the derivates of the error \f$E\f$ with respect to the
							   *  parameters #w and stores the result in the elements of the array
							   *  #dedw. As a default implementation, the derivatives are estimated
 *  numerically by using small variations #epsilon:
							   *  
							   *  \f[\frac{\partial E}{\partial \omega_i} = 
							   *      \frac{E(\omega_i + \epsilon) - E(\omega_i)}{\epsilon} 
							   *      + O(\epsilon) \f]
							   *
							   *  However, to improve the speed, it is adviserable to overload this
							   *  method by a more sophisticated calculation, depending on the
							   *  method #error. Usually, this method is overloaded by members of
							   *  classes that realize certain error measures. Nevertheless, this
							   *  default implementation can be used to check new implementations of
							   *  the derivative of the error with respect to the weights.
							   *  
							   *  Usually, as a byproduct of the calculation of the derivative one gets the 
							   *  the error \f$E\f$ itself very efficiently. Therefore, the method #error
							   *  gives back this value. This additional effect can be switched of by means
							   *  of the third parameter (returnError = false).
							   *  
							   *  
							   *      \param  input Vector of input values.
							   *      \param  target Vector of output values.
							   *      \param  returnError Determines whether or not to calculate the error 
							   *                          itself. 
							   *      \return The error \em E if \em returnError is set to "true", "-1"
							   *              otherwise.
							   *
							   *
							   *  \author  C. Igel
							   *  \date    1999
							   *
							   *  \par Changes
							   *      none
							   *
							   *  \par Status
							   *      stable
							   *
 */
									  virtual double   derror  (const Array<double> &input, const Array<double> &target,
									  bool returnError = true) {
										  double e, eEpsilon, wMem;
										  e = error(input, target);
										  for(unsigned i = 0; i < w.nelem(); i++) {
											  wMem     = w(i);
											  w(i)    += epsilon;
											  eEpsilon = error(input, target);
											  w(i)     = wMem;
											  dedw(i)    = (eEpsilon - e) /  epsilon;
										  }

										  if (!returnError)
											  e=-1;
    
										  return e;  
									  }

//===========================================================================
/*!
									   *  \brief Returns the model's answer #output on the stimulus #input.
												  *
												  *  This method calculates the output of the model depending on the
												  *  #input. The arrays #input and #output can either be one- or
												  *  two-dimensional, depending on whether one or many patterns should
												  *  be processed. The number of elements in the last dimension of
												  *  these arrays must fit the #inputDimension and #outputDimension, in
												  *  case of two-dimensional input, the number of elements in the first
												  *  dimension equals the number of patterns. The method is pure
												  *  virtual as it has to be implemented by the different models (e.g.,
												  *  #FFNet, #RBFNet, #RNNet).
												  *  
												  *      \param  input Vector of input values.
												  *      \param  output Vector of output values.
												  *      \return None.
												  *
												  *  \author  C. Igel
												  *  \date    1999
												  *
												  *  \par Changes
												  *      none
												  *
												  *  \par Status
												  *      stable
												  *
 */
														 virtual void model(const Array<double> &input, Array<double> &output) = 0;

//===========================================================================
/*!
												  *  \brief Calculates the field #dedw, the derivative of the model
															 *  output with respect to the weights.
															 *
															 *  This method performs the calculation of the model output with
															 *  respect to the weights #w. Of course, this calculation can only be
															 *  performed based on a single input pattern; therefore #input must
															 *  be a onedimensional array. As a default implementation, the
															 *  derivatives are estimated numerically by using small variations
 *  #epsilon:
															 *  
															 *  \f[\frac{\partial m_j}{\partial \omega_i} = 
															 *      \frac{m_j(\omega_i + \epsilon) - m_j(\omega_i)}{\epsilon} 
															 *      + O(\epsilon) \f]
															 *
															 *  However, to improve the speed, it is adviserable to overload this
															 *  method by a more sophisticated calculation, depending on the
															 *  method #model. Usually, this method is overloaded by members of
															 *  classes that realize certain models (e.g., FFNet, RBFNet,
															 *  RNNet). Nevertheless, this default implementation can be used to
															 *  check new implementations of the derivative of the model with
															 *  respect to the weights.
															 *  
															 *      \param  input Vector of input values.
															 *      \return None.
															 *
															 *  \author  C. Igel
															 *  \date    1999
															 *
															 *  \par Changes
															 *      none
															 *
															 *  \par Status
															 *      stable
															 *
 */
																	virtual void   dmodel  (const Array<double> &input) {
																double wMem;
																if(input.ndim() == 1) {
																	Array<double> output       (outputDimension);
																	Array<double> outputEpsilon(outputDimension);
																	model(input, output);
																	for(unsigned i = 0; i < w.nelem(); i++) {
																		wMem     = w(i);
																		w(i)    += epsilon;
																		model(input, outputEpsilon);
																		w(i)     = wMem;
																		for(unsigned c = 0; c < outputDimension; c++)
																			dmdw(c, i) = (outputEpsilon(c) - output(c)) / epsilon;
																	}
																} else {
																	std::cerr << "the derivative of the model with respect to more than one input is not defined" << std::endl;
																	exit(EXIT_FAILURE);
																}    
																	}

//===========================================================================
/*!
																	 *  \brief Calculates the field #dedw, the derivative of the model
																				*  output with respect to the weights.
																				*
																				*  This method performs the calculation of the model output with
																				*  respect to the weights #w. Additionally, it calculates the output
																				*  of the model depending on the given input pattern. In many model
																				*  types the model output appears as a byproduct of the calculation
																				*  of the dericative. In this cases it is highly recomanded to
																				*  reimplement #dmodel by such a more efficient implementation.
																				*  
																				*      \param  input Vector of input values.
																				*      \param  output Vector of model output.
																				*      \return None.
																				*
																				*  \author  C. Igel
																				*  \date    1999
																				*
																				*  \par Changes
																				*      none
																				*
																				*  \par Status
																				*      stable
																				*
 */
																				virtual void   dmodel(const Array<double> &input, Array<double> &output) {
																				model (input, output);
																				dmodel(input);
																				}

	protected:
  
  //! Perturbation for gradient approximation. Usually, this value has
  //! to be choosen quite small.
		double epsilon;

//===========================================================================
/*!
		 *  \brief Calculates the model's general derivative.
					*
					*   After a call of #model, #df calculates the gradient \f$\sum_i c_i\,
					*   \frac{\partial}{\partial w^j} f^i_w(x^k)\f$ and stores it in
					*   \em dfdw. The coefficients \f$c_i\f$ can be chosen freely.
					*
					*   The point of this function is that for many models, the complexity
					*   of calculating this linear combination of gradients is the same as
					*   of calculating a single output gradient 
					*   
					*   \f[
					*   \frac{\partial}{\partial \vec w} f^i_w(x^k)
					*   \f]
					*
					*   Via the chain rule, the gradient of _any_
					*   (e.g., error) functional can be calculated.
					*
					*
					*   \em Example \em 1: Let \f$c_i=f^i_w(x)-t^i\f$, then 
					*
					*   \f[
					*   \sum_i c_i\, \frac{\partial}{\partial w^j} f^i_w(x^k) = \frac{\partial}{\partial w^j}\, \Big[ \frac{1}{2}\, \sum_i(f^i_w(x^k)-t^i)^2 \Big]
					*   \f] 
					*
					*   is the parameter gradient of an MSE with respect to target values 
					*   \f$t^i\f$.
					*
					*
					*   \em Example \em 2: Let \f$c_i=\delta_{il}\f$, then the 
 *   gradients of only the \f$l\f$th output will be calculated: 
					*
					*   \f[
					*   \frac{\partial}{\partial w^j} f^l_w(x^k)
					*   \f]
					*
					*   and
					*   
					*   \f[
					*   \frac{\partial}{\partial x^k} f^l_w(x^k)
					*   \f]
					* 
					*   This example is invoked when #dmodel calculates
					*   the gradient of each output with respect to the parameters.
					*
					*   The default implementation is the approximation
					*   df_app(c, dcf_w, dcf_x). 
					*  
					*      \param  c the corefficients \f$c_i\f$ in the formula above.
					*      \param  dfdw the calculated model's general derivative.
					*      \return None.
					*
					*  \author  M. Toussaint, C. Igel
					*  \date    2002
					*
					*  \par Changes
					*      none
					*
					*  \par Status
					*      under construction
					*
 */
						   virtual void df(const Array<double> &c, Array<double> &dfdw) 
				   { 
					   std::cerr << "function df(const Array<double> &c, Array<double> &dfdw) "
							   << "used but not implemented" << std::endl;
					   exit(EXIT_FAILURE);
				   }

};

#endif










// ***************** IT FOLLOWS THE REPLACEMENT VERSION
// ***************** FOR THE NEW RECLAM STRUCTURE
// ***************** PROVIDED FOR DOWNWARDS COMPATIBILITY
#else




//===========================================================================
/*!
 *  \file ModelInterface.h
 *
 *  \brief [DECREPATED] Realizes the communication between the different modules.
 *
 *  <b>ATTENTION: THIS FILE IS DECREPATED!</b><br>
 *  It is provided for downward compatibility with Shark versions
 *  up to 1.4.x only. Programmers are discouraged from using it.
 *  A new design is defined through the classes #Model,
 *  #ErrorFunction and #Optimizer.
 *  This class is a workaround which defines some dummy functions
 *  and data members that mimic the behaviour of the original
 *  ModelInterface class. It's description follows:<br>
 *
 *  To provide flexibility ReClaM offers different modules that can
 *  be put together to form an environment for solving regression
 *  and classification tasks.<BR>
 *  The choice of three module types is necessary: A parametric model,
 *  an error function module and an optimization algorithm model.
 *  All modules exchange parameters and other information, so the task
 *  of this class is to offer this communication.
 *
 *  \author  C. Igel
 *  \date    1999
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
 *      $RCSfile: ModelInterface.h,v $<BR>
 *      $Revision: 2.2 $<BR>
 *      $Date: 2005/12/30 10:50:17 $
 *
 *  \par Changes:
 *      $Log: ModelInterface.h,v $
 *      Revision 2.2  2005/12/30 10:50:17  glasmtbl
 *      *** empty log message ***
 *
 *      Revision 2.0  2003/11/28 16:23:22  shark-admin
 *      Revision tag reset to revision tag 2.x
 *
 *      Revision 1.1.1.1  2003/11/24 13:37:06  shark-admin
 *      INI Administration
 *
 *      Revision 1.6  2002/09/24 08:27:03  rudi
 *      doxygen comments adapted to the new version; new method "df" added.
 *
 *      Revision 1.5  2002/05/16 13:26:03  rudi
 *      doxygen commands added/modified.
 *
 *      Revision 1.4  2001/11/30 14:58:09  rudi
 *      de -> dedw
 *      dw -> dmdw
 *      w, dedw, dmdw protected
 *      fderror removed
 *      derror returns error, third parameter
 *      fdmodel(.,.) -> dmodel(.,.)
 *      regularization removed
 *      doxygen comments added.
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
 *
 */
//===========================================================================

#ifndef _ModelInterface_H_
#define _ModelInterface_H_


#include "Model.h"
#include "ErrorFunction.h"


/*!
 *  <b>ATTENTION: THIS CLASS IS DECREPATED!</b><br>
 *  It is provided for downward compatibility with Shark versions
 *  up to 1.4.x only. Programmers are discouraged from using it.
 *  A new design is defined through the classes #Model,
 *  #ErrorFunction and #Optimizer.
 *  This class is a workaround which defines some dummy functions
 *  and data members that mimic the behaviour of the original
 *  ModelInterface class.<br>
*/
class ModelInterface : public ErrorFunction, public Model
{
public:
	//! The constructor defines references to those variables
	//! whose names have changed in the context of the ReClaM
	//! re-design replacing #ModelInterface by #Model.
	ModelInterface()
		: w(Model::parameter)
	{
	}

	//! Descructor
	virtual ~ModelInterface()
	{
	}


	//! Model output - this call is passed through to #Model.
	virtual void dmodel(const Array<double>& input)
	{
		((Model*)this)->modelDerivative(input, dmdw);
	}

	//! Derivative of the model output - this call is passed through to #Model.
	virtual void dmodel(const Array<double>& input, Array<double>& output)
	{
		((Model*)this)->modelDerivative(input, output, dmdw);
	}

	//! Error evaluation - this call is passed through to #ErrorFunction.
	virtual double error(const Array<double>& input, const Array<double>& target)
	{
		return ((ErrorFunction*)this)->error(*(Model*)this, input, target);
	}

	//! Derivative of the error - this call is passed through to #ErrorFunction.
	virtual double derror(const Array<double>& input, const Array<double>& target, bool bReturnError = true)
	{
		return ((ErrorFunction*)this)->errorDerivative(*(Model*)this, input, target, dedw);
	}

protected:
	//! Reference to #parameter
	Array<double>& w;

	//! Derivative of the model output w.r.t #w, computed by #dmodel
	Array<double> dmdw;

	//! Derivative of the error w.r.t. #w, computed by #derror
	Array<double> dedw;
};


#endif


#endif
