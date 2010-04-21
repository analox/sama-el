//===========================================================================
/*!
*  \file optCMA.h
*
*  \brief The CMA-ES as a ReClaM Optimizer
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
*      $RCSfile: optCMA.h,v $<BR>
*
*  \par Changes:
*      $Log: optCMA.h,v $
*      Revision 1.3  2006/04/10 12:59:49  glasmtbl
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

#ifndef _optCMA_H_
#define _optCMA_H_


#include <ReClaM/Optimizer.h>
#include <EALib/cma2004.h>


class optCMA : public Optimizer
{
public:
	//! Constructor
	optCMA()
	{
		parents = NULL;
		offspring = NULL;
	}

	//! Destructor
	~optCMA()
	{
		if (parents != NULL) delete parents;
		if (offspring != NULL) delete offspring;
	}


	//! basic initialization with default parameters
	//! \param  model  model to optimize
	void init(Model& model)
	{
		init(model, 0.01);
	}

	//! initialization the additional parameter
	//! \param  model  model to optimize
	//! \param  sigma  initial step size
	void init(Model& model, double sigma)
	{
		if (parents != NULL) delete parents;
		if (offspring != NULL) delete offspring;

		int i;
		int dim = model.getParameterDimension();
		int lambda = cma.suggestLambda(dim);
		int mu = cma.suggestMu(lambda);

		ChromosomeT<double> chrom_w(dim);
		ChromosomeT<double> chrom_0(dim);
		for (i=0; i<dim; i++)
		{
			chrom_w[i] = model.getParameter(i);
			chrom_0[i] = 0.0;
		}

		parents = new Population(mu, chrom_w, chrom_0);
		parents->setMinimize();
		offspring = new Population(lambda, chrom_w, chrom_0);
		offspring->setMinimize();

		cma.init(dim, sigma, *parents);
	}

	//! create and select one CMA-ES generation
	double optimize(Model& model, ErrorFunction& errorfunction, const Array<double>& input, const Array<double>& target)
	{
		int i, o;
		double err;
		Individual* pI;
		ChromosomeT<double>* pC;
		int dim = model.getParameterDimension();
		int lambda = offspring->size();

		// create lambda feasible offspring
		for (o=0; o<lambda; o++)
		{
			while (true)
			{
				pI = &(offspring->operator [] (o));
				cma.create(*pI);
				pC = (ChromosomeT<double>*)(&(*pI)[0]);
				for (i=0; i<dim; i++)
				{
					model.setParameter(i, (*pC)[i]);
				}
				if (! model.isFeasible()) continue;
				err = errorfunction.error(model, input, target);
				pI->setFitness(err);
				break;
			}
		}

		parents->selectMuLambda(*offspring, 0);
		cma.updateStrategyParameters(*parents);

		// Copy the best solution to the model and
		// return the corresponding error value.
		pI = &(parents->operator [] (0));
		pC = (ChromosomeT<double>*)(&(*pI)[0]);
		for (i=0; i<dim; i++)
		{
			model.setParameter(i, (*pC)[i]);
		}
		return pI->fitnessValue();
	}

protected:
	//! CMA object from EALib
	CMA2004 cma;

	//! parent population
	Population* parents;

	//! offspring population
	Population* offspring;
};


#endif
