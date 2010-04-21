//===========================================================================
/*!
*  \file ConcatenatedModel.h
*
*  \brief The ConcatenatedModel encapsulates a chain of basic models.
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
*      $RCSfile: ConcatenatedModel.h,v $<BR>
*
*  \par Changes:
*      $Log: ConcatenatedModel.h,v $
*      Revision 1.6  2006/04/20 16:16:51  glasmtbl
*      isFeasible() added
*
*      Revision 1.4  2006/01/16 15:36:30  glasmtbl
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

#ifndef _ConcatenatedModel_H_
#define _ConcatenatedModel_H_


#include <ReClaM/Model.h>
#include <vector>


class ConcatenatedModel : public Model
{
public:
	//! Constructor
	ConcatenatedModel()
	{
	}

	//! Destructor
	~ConcatenatedModel()
	{
		int i, ic = models.size();
		for (i=0; i<ic; i++) delete models[i];
		models.clear();
	}


	//! Append a model to the chain.
	//! Note that a once added model becomes part of the ComcatenatedModel
	//! object, that is the ComcatenatedModel object destroys the added
	//! model in its destructor.
	void AppendModel(Model* pModel)
	{
		if (models.size() > 0)
		{
			Model* pLast = models[models.size() - 1];
			if (pLast->getOutputDimension() != pModel->getInputDimension())
			{
				throw("[ConcatenatedModel::AppendModel] Dimension conflict");
			}

			int d = parameter.dim(0);
			int i, ic = pModel->getParameterDimension();
			parameter.resize(d + ic, true);
			for (i=0; i<ic; i++) parameter(d + i) = pModel->getParameter(i);
		}
		else
		{
			int i, ic = pModel->getParameterDimension();
			parameter.resize(ic, false);
			for (i=0; i<ic; i++) parameter(i) = pModel->getParameter(i);
			inputDimension = pModel->getInputDimension();
		}
		outputDimension = pModel->getOutputDimension();

		models.push_back(pModel);
	}

	//! Evaluate the concatenated model
	void model(const Array<double>& input, Array<double> &output)
	{
		Array<double> temp[2];
		int i, ic = models.size();
		int t = 0;
		int p, pc, pp = 0;

		// copy the parameters into the elementary models
		for (i=0; i<ic; i++)
		{
			pc = models[i]->getParameterDimension();
			for (p=0; p<pc; p++)
			{
				models[i]->setParameter(p, parameter(pp));
				pp++;
			}
		}

		// propagate the input through the model chain
		for (i=0; i<ic; i++)
		{
			if (i == 0)
			{
				models[i]->model(input, temp[0]);
			}
			else if (i == ic - 1)
			{
				models[i]->model(temp[t], output);
			}
			else
			{
				models[i]->model(temp[t], temp[1 - t]);
				t = 1 - t;
			}
		}
	}

	//! We will have to implement the chain rule for all cases...
	void modelDerivative(const Array<double>& input, Array<double>& derivative)
	{
		int i, ic = models.size();
		int t = 0;
		int p, pc, pp = 0;

		// copy the parameters into the elementary models
		for (i=0; i<ic; i++)
		{
			pc = models[i]->getParameterDimension();
			for (p=0; p<pc; p++)
			{
				models[i]->setParameter(p, parameter(pp));
				pp++;
			}
		}

		// TODO: this means some work not done right now ...
	}

	//! We will have to implement the chain rule for all cases...
	void modelDerivative(const Array<double>& input, Array<double>& output, Array<double>& derivative)
	{
		int i, ic = models.size();
		int t = 0;
		int p, pc, pp = 0;

		// copy the parameters into the elementary models
		for (i=0; i<ic; i++)
		{
			pc = models[i]->getParameterDimension();
			for (p=0; p<pc; p++)
			{
				models[i]->setParameter(p, parameter(pp));
				pp++;
			}
		}

		// TODO: perform both model and modelDerivative at the same time ...
	}

	//! Check whether the parameters define a valid model
	bool isFeasible()
	{
		int i, ic = models.size();
		if (! models[i]->isFeasible()) return false;
		return true;
	}

protected:
	//! chain of models
	std::vector<Model*> models;
};


#endif
