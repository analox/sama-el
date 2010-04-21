//===========================================================================
/*!
*  \file LinearModel.h
*
*  \brief Linear function on a real vector space
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
*      $RCSfile: LinearModel.h,v $<BR>
*
*  \par Changes:
*      $Log: LinearModel.h,v $
*      Revision 1.3  2006/01/27 09:49:40  glasmtbl
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

#ifndef _LinearModel_H_
#define _LinearModel_H_


#include <ReClaM/Model.h>


class LinearModel : public Model
{
public:
	LinearModel(int dimension)
	{
		inputDimension = dimension;
		outputDimension = 1;

		parameter.resize(dimension);
		parameter = 0.0;
	}

	~LinearModel()
	{
	}


	void model(const Array<double>& input, Array<double> &output)
	{
		if (input.ndim() == 1)
		{
			int i;
			double value;
			output.resize(1, false);
			value = 0.0;
			for (i=0; i<inputDimension; i++) value += input(i) * parameter(i);
			output(0) = value;
		}
		else if (input.ndim() == 2)
		{
			int j, jc = input.dim(0);
			int i;
			double value;
			output.resize(jc, false);
			for (j=0; j<jc; j++)
			{
				value = 0.0;
				for (i=0; i<inputDimension; i++) value += input(j, i) * parameter(i);
				output(j) = value;
			}
		}
		else throw("[LinearModel::model] invalid number of dimensions.");
	}

	void modelDerivative(const Array<double>& input, Array<double>& derivative)
	{
		if (input.ndim() == 1)
		{
			derivative = parameter;
		}
		else if (input.ndim() == 2)
		{
			int j, jc = input.dim(0);
			derivative.resize(jc, getParameterDimension());
			for (j=0; j<jc; j++) derivative[j] = parameter;
		}
		else throw("[LinearModel::modelDerivative] invalid number of dimensions.");
	}
};


#endif
