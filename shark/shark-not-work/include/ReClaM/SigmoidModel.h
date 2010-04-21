//===========================================================================
/*!
*  \file SigmoidModel.h
*
*  \brief sigmoid map \f$ x \mapsto \frac{1}{1+exp(Ax+B)} \f$
*
*  \author  T. Glasmachers
*  \date    2006
*
*  \par Copyright (c) 1999-2006:
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
*      $RCSfile: SigmoidModel.h,v $<BR>
*
*  \par Changes:
*      $Log: SigmoidModel.h,v $
*      Revision 2.1  2006/01/27 09:49:40  glasmtbl
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

#ifndef _SigmoidModel_H_
#define _SigmoidModel_H_


#include <ReClaM/Model.h>
#include <ReClaM/SigmoidModel.h>


class SigmoidModel : public Model
{
public:
	SigmoidModel(double A = 1.0, double B = 0.0)
	{
		inputDimension = 1;
		outputDimension = 1;

		parameter.resize(2);
		parameter(0) = A;
		parameter(1) = B;
	}

	~SigmoidModel()
	{
	}


	inline double get_A()
	{
		return parameter(0);
	}

	inline double get_B()
	{
		return parameter(1);
	}

	void model(const Array<double>& input, Array<double> &output)
	{
		if (input.ndim() == 1)
		{
			SIZE_CHECK(input.dim(0) == 1);
			output.resize(1, false);
			output(0) = 1.0 / (1.0 + exp(parameter(0) * input(0) + parameter(1)));
		}
		else if (input.ndim() == 2)
		{
			int j, jc = input.dim(0);
			SIZE_CHECK(input.dim(1) == 1);
			output.resize(jc, false);
			for (j=0; j<jc; j++)
			{
				output(j) = 1.0 / (1.0 + exp(parameter(0) * input(j, 0) + parameter(1)));
			}
		}
		else throw("[SigmoidModel::model] invalid number of dimensions.");
	}

	void modelDerivative(const Array<double>& input, Array<double>& derivative)
	{
		if (input.ndim() == 1)
		{
			SIZE_CHECK(input.dim(0) == 1);
			derivative.resize(2, false);
			double e = exp(parameter(0) * input(0) + parameter(1));
			double f = 1.0 / (1.0 + e);
			double d = -f*f*e;
			derivative(0) = d*input(0);
			derivative(1) = d;
		}
		else if (input.ndim() == 2)
		{
			int j, jc = input.dim(0);
			SIZE_CHECK(input.dim(1) == 1);
			derivative.resize(jc, 2, false);
			for (j=0; j<jc; j++)
			{
				double e = exp(parameter(0) * input(0) + parameter(1));
				double f = 1.0 / (1.0 + e);
				double d = -f*f*e;
				derivative(j, 0) = d*input(0);
				derivative(j, 1) = d;
			}
		}
		else throw("[SigmoidModel::modelDerivative] invalid number of dimensions.");
	}

	void modelDerivative(const Array<double>& input, Array<double>& output, Array<double>& derivative)
	{
		if (input.ndim() == 1)
		{
			SIZE_CHECK(input.dim(0) == 1);
			derivative.resize(2, false);
			output.resize(1, false);
			double e = exp(parameter(0) * input(0) + parameter(1));
			double f = 1.0 / (1.0 + e);
			double d = -f*f*e;
			output(0) = f;
			derivative(0) = d*input(0);
			derivative(1) = d;
		}
		else if (input.ndim() == 2)
		{
			int j, jc = input.dim(0);
			SIZE_CHECK(input.dim(1) == 1);
			derivative.resize(jc, 2, false);
			output.resize(jc, false);
			for (j=0; j<jc; j++)
			{
				double e = exp(parameter(0) * input(0) + parameter(1));
				double f = 1.0 / (1.0 + e);
				double d = -f*f*e;
				output(j) = f;
				derivative(j, 0) = d*input(0);
				derivative(j, 1) = d;
			}
		}
		else throw("[SigmoidModel::modelDerivative] invalid number of dimensions.");
	}
};


#endif
