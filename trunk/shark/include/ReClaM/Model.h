//===========================================================================
/*!
*  \file Model.h
*
*  \brief Base class of all models.
*
* ReClaM provides the three base classes Model, ErrorFunction and
* Optimizer which make up the ReClaM framework for solving regression
* and classification tasks. This design overrides the ModelInterface
* design which is kept for downward compatibility.<BR>
* The Model class encapsulates a data processing trainable model.
* Its internal state is described by an array of parameters.
*
*  \author  T. Glasmachers
*  \date    2005
*
*  \par Copyright (c) 1999-2005:
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
*      $RCSfile: Model.h,v $<BR>
*
*  \par Changes:
*      $Log: Model.h,v $
*      Revision 2.8  2006/04/20 11:34:14  glasmtbl
*      *** empty log message ***
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
*
*/
//===========================================================================

#ifndef _Model_H_
#define _Model_H_


#include <Array/Array.h>
#include <fstream>


// forward declarations
class Optimizer;
class ErrorFunction;


class Model
{
public:
	//! Constructor
	Model()
	{
		epsilon = 1e-8;
	}

	//! Destructor
	virtual ~Model()
	{
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
	*  \author  T. Glasmachers
	*  \date    2005
	*
	*  \par Changes
	*      none
	*
	*  \par Status
	*      stable
	*
	*/
	virtual void model(const Array<double>& input, Array<double>& output) = 0;

	//! Do a model evaluation on a const object
	void model(const Array<double>& input, Array<double> &output) const
	{
		Model* pT = const_cast<Model*>(this);
		return pT->model(input, output);
	}

	//===========================================================================
	/*!
	*  \brief Calculates the field #dmdw, the derivative of the model
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
	*      \param  derivative Matrix of partial derivatives.
	*      \return None.
	*
	*  \author  T. Glasmachers
	*  \date    2005
	*
	*  \par Changes
	*      none
	*
	*  \par Status
	*      stable
	*
	*/
	virtual void modelDerivative(const Array<double>& input, Array<double>& derivative)
	{
		double old;
		int p, pc = parameter.dim(0);
		Array<double> output;
		Array<double> perturbed_output;
		model(input, output);
		int o, oc = output.dim(0);
		derivative.resize(oc, pc, false);
		for (p=0; p<pc; p++)
		{
			old = parameter(p);
			parameter(p) += epsilon;
			model(input, perturbed_output);
			for (o=0; o<oc; o++)
			{
				derivative(o, p) = (perturbed_output(o) - output(o)) / epsilon;
			}
			parameter(p) = old;
		}
	}

	//! Compute the model derivative on a const object
	void modelDerivative(const Array<double>& input, Array<double>& derivative) const
	{
		Model* pT = const_cast<Model*>(this);
		return pT->modelDerivative(input, derivative);
	}

	//===========================================================================
	/*!
	*  \brief Calculates the field derivative, the derivative of the model
	*  output with respect to the parameters.
	*
	*  This method performs the calculation of the model output with
	*  respect to the parameters #parameter. Additionally, it calculates the
	*  output of the model depending on the given input pattern. In many model
	*  types the model output appears as a byproduct of the calculation
	*  of the dericative. In these cases it is highly recommended to
	*  reimplement #dmodel by such a more efficient implementation.
	*
	*      \param  input Vector of input values.
	*      \param  output Vector of model output.
	*      \param  derivative Matrix of partial derivatives.
	*      \return None.
	*
	*  \author  T. Glasmachers
	*  \date    2005
	*
	*  \par Changes
	*      none
	*
	*  \par Status
	*      stable
	*
	*/
	virtual void modelDerivative(const Array<double>& input, Array<double>& output, Array<double>& derivative)
	{
		model(input, output);
		modelDerivative(input, derivative);
	}

	//! Compute the model derivative on a const object
	void modelDerivative(const Array<double>& input, Array<double>& output, Array<double>& derivative) const
	{
		Model* pT = const_cast<Model*>(this);
		return pT->modelDerivative(input, output, derivative);
	}

	//===========================================================================
	/*!
	*  \brief Calculates the model's general derivative.
	*
	*   After a call of #model, #generalDerivative calculates the gradient \f$\sum_i c_i\,
	*   \frac{\partial}{\partial w^j} f^i_w(x^k)\f$ and stores it in
	*   \em derivative. The coefficients \f$c_i\f$ can be chosen freely.
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
	*   This example is invoked when #modelDerivative calculates
	*   the gradient of each output with respect to the parameters.
	*
	*   The default implementation is the computationally inefficient
	*   usage of modelDerivative.
	*
	*      \param  input the corefficients \f$c_i\f$ in the formula above.
	*      \param  coefficient the corefficients \f$c_i\f$ in the formula above.
	*      \param  derivative the calculated model's general derivative.
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
	virtual void generalDerivative(const Array<double>& input, const Array<double>& coefficient, Array<double>& derivative)
	{
		Array<double> md;
		modelDerivative(input, md);

		int i, ic = getOutputDimension();
		int j, jc = getParameterDimension();

		derivative.resize(jc, false);

		for (j=0; j<jc; j++)
		{
			double td = 0.0;
			for (i=0; i<ic; i++)
			{
				td += coefficient(i) * md(i, j);
			}
			derivative(j) = td;
		}
	}

	/*!
	*  \brief check whether the parameters define a feasible model
	*
	*  The default implementation returns true, that is,
	*  every parameter configuration is considered feasible
	*  and unconstrained optimization is applicable.
	*  It is the Optimizer's responsibility to check the
	*  isFeasible() flag.
	*  
	*     \return true if the model is feasible, false otherwise
	*
	*  \author  T. Glasmachers
	*  \date    2006
	*
	*  \par Changes
	*      none
	*
	*  \par Status
	*      stable
	*
	*/
	virtual bool isFeasible()
	{
		return true;
	}

	//! check whether the parameters define a feasible model on a const object
	bool isFeasible() const
	{
		Model* pT = const_cast<Model*>(this);
		return pT->isFeasible();
	}

	//! Returns the dimension of the model input.
	const inline unsigned int getInputDimension() const
	{
		return inputDimension;
	}

	//! Returns the dimension of the model output.
	const inline unsigned int getOutputDimension() const
	{
		return outputDimension;
	}

	//! Returns the number of optimizable model parameters, i.e. the dimension of parameter.
	const inline unsigned int getParameterDimension() const
	{
		SIZE_CHECK(parameter.ndim() == 1);
		return parameter.dim(0);
	}

	//! Returns a specific model parameter.
	virtual double getParameter(unsigned int index) const
	{
		return parameter(index);
	}

	//! Modifies a specific model parameter.
	virtual void setParameter(unsigned int index, double value)
	{
		parameter(index) = value;
	}

	//===========================================================================
	/*!
	 *  \brief Sets the value of #epsilon.
	 *
	 *  The value of #epsilon is utilized for a numeric estimation of the
	 *  derivative of the error, with respect to the model parameters.
	 *  #epsilon should be choosen as a small value. The default value
	 *  is #epsilon = 1e-8.
	 *
	 *      \param  epsilon New value for epsilon.
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
	inline void setEpsilon (double eps)
	{ 
		epsilon = eps;
	};

protected:
	//! Model parameter vector.
	//! #parameter is a 1-dimensional array with the number
	//! of elements equal to the number of free parameters of
	//! the model (i.e., equal to the number of weights of a
	//! neural network).
	Array<double> parameter;

	//! Dimension of the data accepted by the model as input.
	unsigned int inputDimension;

	//! Dimension of the model output per single input pattern.
	unsigned int outputDimension;

	//! Precision parameter for the numerical error gradient approximation.
	double epsilon;

public:
	//! \brief Read the model parameters from a stream.
	//!
	//! The model parameters are read from a text stream.
	//! The single numbers are separated with a single
	//! space character. The line break CR/LF is used as
	//! an end marker. The setParameter member is used
	//! to set the parameter values in order to support
	//! models that override the parameter manegement.
	//!
	//! \author  T. Glasmachers
	//! \date    2006
	//!
	friend inline std::istream& operator >> (std::istream& is, Model& model)
	{
		char c;
		double value;
		std::vector<double> v;
		int p, pc;
		while (true)
		{
			is >> value;
			v.push_back(value);
			is.read(&c, 1);
			if (c == '\r')
			{
				is.read(&c, 1);
				if (! is.good() || c != '\n') throw "[Model::operator >>] invalid data format";
				break;
			}
			if (! is.good() || c != ' ') throw "[Model::operator >>] invalid data format";
		}
		pc = v.size();
		model.parameter.resize(pc, false);
		for (p=0; p<pc; p++) model.setParameter(p, v[p]);
		return is;
	}

	//! \brief Write the model parameters to a stream.
	//!
	//! The model parameters are written to a text stream.
	//! The single numbers are separated with a single
	//! space character. The line break CR/LF is used as
	//! an end marker.
	//!
	//! \author  T. Glasmachers
	//! \date    2006
	//!
	friend inline std::ostream& operator << (std::ostream& os, const Model& model)
	{
		int p, pc = model.parameter.dim(0);
		char buffer[50];
		for (p=0; p<pc; p++)
		{
			if (p != 0) os.write(" ", 1);
			sprintf(buffer, "%.15g", model.parameter(p));
			os.write(buffer, strlen(buffer));
		}
		os.write("\r\n", 2);
		return os;
	}

	//! Read the model parameters from a file.
	//! \author  T. Glasmachers
	//! \date    2006
	bool Load(char* filename)
	{
		std::ifstream is;
		is.open(filename);
		if (! is.is_open()) return false;
		is >> *this;
		is.close();
		return true;
	}

	//! Write the model parameters to a file.
	//! \author  T. Glasmachers
	//! \date    2006
	bool Save(char* filename)
	{
		std::ofstream os;
		os.open(filename);
		if (! os.is_open()) return false;
		os << *this;
		os.close();
		return true;
	}
};


#endif
